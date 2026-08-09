#!/usr/bin/env julia

# ==============================================================================
# Dynamic mean-field theory (DMFT) solver for a random rate network with a
# nonlinear self-interaction
# ==============================================================================
#
# Microscopic network:
#
#     dx_i/dt = -x_i + s*tanh(x_i)
#                     + g * sum_{j != i} J_ij*tanh(x_j),
#
# where J_ij are independent random variables with
#
#     E[J_ij] = 0,       E[J_ij^2] = 1/N.
#
# For fully asymmetric, independent couplings, the N -> infinity DMFT is the
# single-site stochastic equation
#
#     dx/dt = -x + s*tanh(x) + eta(t),                         (DMFT-1)
#
# where eta(t) is a zero-mean Gaussian process whose covariance must satisfy
#
#     E[eta(t) eta(t')] = g^2 C(t-t'),                         (DMFT-2)
#
#     C(t-t') = E[tanh(x(t)) tanh(x(t'))].                     (DMFT-3)
#
# Thus C is a functional fixed point.  Given a trial C, we draw Gaussian
# colored-noise paths with covariance g^2 C, integrate (DMFT-1), estimate the
# right side of (DMFT-3), and repeat until the input and output correlations
# agree.
#
# The numerical method below is the Fourier/Monte-Carlo iteration described in
# the paper on strongly self-coupled random networks:
#
#   1. Represent a stationary covariance C(tau) on a periodic time grid.
#   2. Diagonalize its circulant covariance matrix with an FFT.
#   3. Draw many Gaussian eta(t) paths with covariance g^2 C(tau).
#   4. Integrate the nonlinear scalar equation for every path.
#   5. Average tanh(x(t))*tanh(x(t+tau)) over paths and time origins.
#   6. Mix this new estimate with the previous C and iterate.
#
# Important scope condition:
# This code assumes independent, fully asymmetric J_ij.  If J_ij and J_ji are
# correlated, the effective process also contains a retarded response/memory
# term, and the equations solved here are no longer complete.
#
# Dependencies:
#   * Julia standard libraries: Random, LinearAlgebra, Printf
#   * FFTW.jl
#
# Install FFTW once from the Julia REPL with
#
#     import Pkg
#     Pkg.add("FFTW")
#
# Run this file with
#
#     julia nonlinear_self_dmft.jl
#
# The example at the bottom writes CSV files containing the converged
# correlation, convergence history, spectrum, and a few sample paths.
# ============================================================================== 

module NonlinearSelfDMFT

using FFTW
using LinearAlgebra
using Printf
using Random

export DMFTParams,
       DMFTResult,
       solve_dmft,
       save_result,
       positive_well_location


# ------------------------------------------------------------------------------
# User-adjustable numerical parameters
# ------------------------------------------------------------------------------

Base.@kwdef struct DMFTParams
    # Model parameters in
    #     dx/dt = -x + s*tanh(x) + eta(t)
    # with Cov[eta] = g^2 C.
    s::Float64 = 1.2
    g::Float64 = 1.5

    # Time discretization.  The represented period is T = n_time*dt.
    # A small dt controls ODE integration error; a large T is needed whenever
    # C(tau) decays slowly, especially near a transition or at large s.
    dt::Float64 = 0.05
    n_time::Int = 4096

    # Number of independent effective single-site trajectories used in the
    # Monte-Carlo average in (DMFT-3).
    n_paths::Int = 128

    # The sampled colored noise is periodic on the numerical time window.
    # Before recording a path, we integrate through this many whole periods so
    # that the arbitrary initial value of x is largely forgotten.  In a
    # bistable regime, genuine dependence on the initial basin can remain.
    warmup_periods::Int = 2

    # Fixed-point iteration settings:
    #
    #   C_next = (1-mixing)*C_old + mixing*C_output.
    #
    # Smaller mixing is slower but more robust when the DMFT map is stiff.
    mixing::Float64 = 0.20
    max_iterations::Int = 100
    min_iterations::Int = 10

    # Convergence is declared when either the relative RMS residual or the
    # absolute RMS residual is below its tolerance.  The absolute criterion is
    # important when the solution is C approximately 0.
    relative_tolerance::Float64 = 2.0e-3
    absolute_tolerance::Float64 = 5.0e-4

    # Initial trial covariance:
    #
    #   C_initial(tau) = amplitude * exp(-d_periodic(tau)/time_scale).
    #
    # It is subsequently projected onto the cone of valid positive-semidefinite
    # circulant covariance matrices.
    initial_covariance_amplitude::Float64 = 0.50
    initial_covariance_time_scale::Float64 = 5.0

    # Initial state used for every scalar ODE path.
    # Available choices:
    #   :auto      -> :bimodal for s>1, otherwise :gaussian
    #   :zero      -> x(0)=0
    #   :gaussian  -> x(0)~Normal(0, initial_x_scale^2)
    #   :bimodal   -> equal numbers near +/-x_star, where x_star=s*tanh(x_star)
    #   :positive  -> all paths start near +x_star
    #   :negative  -> all paths start near -x_star
    #
    # In a bistable region, different choices can select different DMFT
    # solutions.  This is not merely a numerical nuisance: static and dynamic
    # solutions can coexist.
    initial_x_mode::Symbol = :auto
    initial_x_scale::Float64 = 1.0

    # Reusing the same underlying white-noise samples at every DMFT iteration
    # makes the finite-M Monte-Carlo map deterministic and usually improves
    # convergence.  Set false to redraw the Monte-Carlo ensemble each iteration.
    common_random_numbers::Bool = true
    seed::Int = 12345

    # Number of eta(t), x(t), and tanh(x(t)) traces retained in the result.
    n_store_paths::Int = 4

    verbose::Bool = true
end


# ------------------------------------------------------------------------------
# Returned data
# ------------------------------------------------------------------------------

struct DMFTResult
    params::DMFTParams

    # Full periodic covariance arrays.  Index 1 is lag zero; index k+1 is lag
    # k*dt.  Because the representation is periodic, C[N-k] is the negative-lag
    # counterpart of C[k].
    C::Vector{Float64}
    C_from_paths::Vector{Float64}
    noise_covariance::Vector{Float64}

    # Nonnegative lags 0, dt, ..., T/2.
    lags::Vector{Float64}

    # FFT-order angular frequencies and continuous-time-style spectral
    # densities dt*FFT(C).  The FFT eigenvalues used internally do not include
    # the factor dt; this output does, so it approximates an integral Fourier
    # transform as dt -> 0.
    angular_frequencies::Vector{Float64}
    C_spectrum::Vector{Float64}
    noise_spectrum::Vector{Float64}

    # Convergence diagnostics, one entry per fixed-point iteration.
    relative_residual_history::Vector{Float64}
    absolute_residual_history::Vector{Float64}
    step_history::Vector{Float64}
    C0_history::Vector{Float64}
    mean_rate_history::Vector{Float64}

    converged::Bool
    iterations::Int
    final_relative_residual::Float64
    final_absolute_residual::Float64
    final_mean_rate::Float64

    # A few final effective-process paths.  Each matrix has size
    # n_time x n_store_paths.
    sample_eta::Matrix{Float64}
    sample_x::Matrix{Float64}
    sample_rate::Matrix{Float64}
end


# ------------------------------------------------------------------------------
# Basic validation
# ------------------------------------------------------------------------------

function validate_params(p::DMFTParams)
    p.dt > 0.0 || throw(ArgumentError("dt must be positive"))
    p.n_time >= 8 || throw(ArgumentError("n_time must be at least 8"))
    iseven(p.n_time) || throw(ArgumentError(
        "n_time must be even in this implementation so that T/2 is a grid lag"
    ))
    p.n_paths >= 1 || throw(ArgumentError("n_paths must be positive"))
    p.warmup_periods >= 0 || throw(ArgumentError("warmup_periods cannot be negative"))
    0.0 < p.mixing <= 1.0 || throw(ArgumentError("mixing must lie in (0,1]"))
    p.max_iterations >= 1 || throw(ArgumentError("max_iterations must be positive"))
    1 <= p.min_iterations <= p.max_iterations || throw(ArgumentError(
        "min_iterations must be between 1 and max_iterations"
    ))
    p.relative_tolerance > 0.0 || throw(ArgumentError(
        "relative_tolerance must be positive"
    ))
    p.absolute_tolerance > 0.0 || throw(ArgumentError(
        "absolute_tolerance must be positive"
    ))
    0.0 <= p.initial_covariance_amplitude <= 1.0 || throw(ArgumentError(
        "initial_covariance_amplitude should lie in [0,1] because |tanh(x)|<=1"
    ))
    p.initial_covariance_time_scale > 0.0 || throw(ArgumentError(
        "initial_covariance_time_scale must be positive"
    ))
    p.initial_x_scale >= 0.0 || throw(ArgumentError(
        "initial_x_scale cannot be negative"
    ))
    0 <= p.n_store_paths <= p.n_paths || throw(ArgumentError(
        "n_store_paths must lie between 0 and n_paths"
    ))

    allowed_modes = (:auto, :zero, :gaussian, :bimodal, :positive, :negative)
    p.initial_x_mode in allowed_modes || throw(ArgumentError(
        "initial_x_mode must be one of $(allowed_modes)"
    ))

    return nothing
end


# ------------------------------------------------------------------------------
# Covariance utilities
# ------------------------------------------------------------------------------

"""
    circular_reverse(c)

Return the covariance sequence with lag tau replaced by -tau on a periodic
array.  If c[1] is lag zero, the returned order is

    [c(0), c(-dt), c(-2dt), ..., c(-(N-1)dt)].
"""
function circular_reverse(c::AbstractVector{<:Real})
    n = length(c)
    out = Vector{Float64}(undef, n)
    out[1] = Float64(c[1])

    @inbounds for k in 2:n
        # k=2 maps to c[n], k=3 maps to c[n-1], ..., k=n maps to c[2].
        out[k] = Float64(c[n - k + 2])
    end

    return out
end


"""
    symmetrize_circular_covariance(c)

Enforce C(tau)=C(-tau), which a real stationary covariance must satisfy.
"""
function symmetrize_circular_covariance(c::AbstractVector{<:Real})
    return 0.5 .* (Float64.(c) .+ circular_reverse(c))
end


"""
    project_to_circulant_covariance(c)

Project a real periodic sequence onto the set of valid circulant covariance
matrices:

  * enforce evenness C(tau)=C(-tau);
  * FFT to obtain the eigenvalues of the circulant covariance matrix;
  * clip negative eigenvalues to zero;
  * transform back.

The clipping is useful for an arbitrary initial guess and removes tiny negative
spectral values caused by roundoff.
"""
function project_to_circulant_covariance(c::AbstractVector{<:Real})
    c_even = symmetrize_circular_covariance(c)

    # For a circulant covariance matrix, FFT(c_even) contains its eigenvalues.
    eigenvalues = real.(fft(c_even))
    eigenvalues .= max.(eigenvalues, 0.0)

    c_projected = real.(ifft(eigenvalues))
    return symmetrize_circular_covariance(c_projected)
end


"""
    initial_covariance(p)

Build a smooth, periodic, exponentially decaying starting guess for C(tau).
"""
function initial_covariance(p::DMFTParams)
    n = p.n_time
    c = Vector{Float64}(undef, n)

    @inbounds for k in 0:(n - 1)
        # On a periodic grid, k and n-k represent opposite directions around
        # the same circle.  Their minimum is the shortest periodic lag.
        periodic_lag = min(k, n - k) * p.dt
        c[k + 1] = p.initial_covariance_amplitude *
                   exp(-periodic_lag / p.initial_covariance_time_scale)
    end

    return project_to_circulant_covariance(c)
end


"""
    circular_autocorrelation(y)

Compute

    A[k] = (1/N) sum_n y[n] y[n+k]

with periodic indexing.  The Wiener-Khinchin identity makes this an O(N log N)
FFT calculation.  No sample mean is subtracted: DMFT requires the full second
moment E[tanh(x(t)) tanh(x(t+tau))], not a connected covariance.
"""
function circular_autocorrelation(y::AbstractVector{<:Real})
    n = length(y)
    yhat = fft(y)
    return real.(ifft(abs2.(yhat))) ./ n
end


# ------------------------------------------------------------------------------
# Initial conditions for the effective scalar process
# ------------------------------------------------------------------------------

"""
    positive_well_location(s)

For s>1, return the positive nonzero root x_star of

    x_star = s*tanh(x_star).

For s<=1 the only root connected to the origin is zero, so return 0.
A bisection method is used to avoid an additional root-finding dependency.
"""
function positive_well_location(s::Real)
    s_float = Float64(s)
    s_float > 1.0 || return 0.0

    f(x) = x - s_float * tanh(x)

    # Just to the right of zero, f(x) is negative because f'(0)=1-s<0.
    lo = sqrt(eps(Float64))

    # Since tanh(x)<1, f(s+1)>1, so this is safely to the right of the root.
    hi = s_float + 1.0

    f(lo) < 0.0 || error("Internal root bracket failure at the lower endpoint")
    f(hi) > 0.0 || error("Internal root bracket failure at the upper endpoint")

    for _ in 1:100
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0
            hi = mid
        else
            lo = mid
        end
    end

    return 0.5 * (lo + hi)
end


"""
    make_initial_states(p, rng)

Generate one initial x value for each effective-process path.
"""
function make_initial_states(p::DMFTParams, rng::AbstractRNG)
    mode = if p.initial_x_mode == :auto
        p.s > 1.0 ? :bimodal : :gaussian
    else
        p.initial_x_mode
    end

    x0 = zeros(Float64, p.n_paths)

    if mode == :zero
        return x0

    elseif mode == :gaussian
        randn!(rng, x0)
        x0 .*= p.initial_x_scale
        return x0

    elseif mode == :bimodal
        x_star = positive_well_location(p.s)

        if x_star == 0.0
            # :bimodal has no distinct wells for s<=1; use a symmetric Gaussian
            # ensemble rather than silently returning an all-zero state.
            randn!(rng, x0)
            x0 .*= p.initial_x_scale
            return x0
        end

        # Pair opposite starting states to represent both basins symmetrically.
        # A small paired jitter avoids putting every path at exactly one value.
        pair_count = div(p.n_paths, 2)
        jitter_scale = 0.05 * p.initial_x_scale

        for pair in 1:pair_count
            jitter = jitter_scale * randn(rng)
            x0[2 * pair - 1] =  x_star + jitter
            x0[2 * pair]     = -x_star - jitter
        end

        if isodd(p.n_paths)
            x0[end] = x_star
        end

        return x0

    elseif mode == :positive || mode == :negative
        x_star = positive_well_location(p.s)
        sign_value = mode == :positive ? 1.0 : -1.0

        # If s<=1, x_star=0.  The optional scale then gives a small displacement
        # in the requested direction.
        base = x_star > 0.0 ? x_star : p.initial_x_scale
        x0 .= sign_value * base
        return x0

    else
        error("Unreachable initial-state mode: $mode")
    end
end


# ------------------------------------------------------------------------------
# Gaussian colored-noise generation
# ------------------------------------------------------------------------------

"""
    make_white_noise_ffts(p, rng)

Precompute FFTs of real, independent, unit-variance white-noise paths.  Reusing
these arrays at every DMFT iteration implements common random numbers.
"""
function make_white_noise_ffts(p::DMFTParams, rng::AbstractRNG)
    white_ffts = Matrix{ComplexF64}(undef, p.n_time, p.n_paths)
    z = Vector{Float64}(undef, p.n_time)

    for path in 1:p.n_paths
        randn!(rng, z)
        white_ffts[:, path] .= fft(z)
    end

    return white_ffts
end


"""
    colored_noise_from_white_fft(sqrt_eigenvalues, white_fft)

Given the square roots of the eigenvalues of a circulant covariance matrix D,
construct a real Gaussian vector eta with Cov(eta)=D:

    eta = IFFT( sqrt(lambda_D) .* FFT(z) ),

where z is real white noise.  With Julia/FFTW's FFT normalization, no additional
factor of N or sqrt(N) is needed.
"""
function colored_noise_from_white_fft(
    sqrt_eigenvalues::AbstractVector{<:Real},
    white_fft::AbstractVector{<:Complex},
)
    return real.(ifft(sqrt_eigenvalues .* white_fft))
end


# ------------------------------------------------------------------------------
# Integrating the nonlinear single-site DMFT equation
# ------------------------------------------------------------------------------

@inline function local_drift(x::Float64, eta::Float64, s::Float64, phi)
    # This is the right side of DMFT-1:
    #
    #     f(x,t) = -x + s*phi(x) + eta(t).
    return -x + s * phi(x) + eta
end


@inline function heun_step(
    x::Float64,
    eta_now::Float64,
    eta_next::Float64,
    dt::Float64,
    s::Float64,
    phi,
)
    # Explicit trapezoidal / Heun step.  The predictor uses eta at t_n and the
    # corrector uses eta at t_{n+1}.  This is second order for a sufficiently
    # smooth sampled forcing and is more accurate than forward Euler at the same
    # dt while remaining transparent.
    f_now = local_drift(x, eta_now, s, phi)
    x_predict = x + dt * f_now
    f_next = local_drift(x_predict, eta_next, s, phi)

    return x + 0.5 * dt * (f_now + f_next)
end


"""
    integrate_effective_path(eta, p, x_initial, phi)

Integrate

    dx/dt = -x + s*phi(x) + eta(t)

for one periodic colored-noise realization.  The same noise period is traversed
`warmup_periods` times; the following period is recorded.
"""
function integrate_effective_path(
    eta::AbstractVector{<:Real},
    p::DMFTParams,
    x_initial::Real,
    phi,
)
    n = p.n_time
    length(eta) == n || throw(DimensionMismatch(
        "eta has length $(length(eta)), but n_time=$n"
    ))

    x = Float64(x_initial)

    # Warm-up passes through the periodic forcing.
    for _ in 1:p.warmup_periods
        @inbounds for k in 1:n
            next_k = (k == n) ? 1 : (k + 1)
            x = heun_step(
                x,
                Float64(eta[k]),
                Float64(eta[next_k]),
                p.dt,
                p.s,
                phi,
            )
        end
    end

    # Record one complete period after warm-up.
    x_trace = Vector{Float64}(undef, n)

    @inbounds for k in 1:n
        x_trace[k] = x
        next_k = (k == n) ? 1 : (k + 1)
        x = heun_step(
            x,
            Float64(eta[k]),
            Float64(eta[next_k]),
            p.dt,
            p.s,
            phi,
        )
    end

    all(isfinite, x_trace) || error(
        "The effective trajectory contains Inf or NaN.  Reduce dt and rerun."
    )

    return x_trace
end


# ------------------------------------------------------------------------------
# One application of the DMFT covariance map C -> F[C]
# ------------------------------------------------------------------------------

"""
    evaluate_dmft_map(C, p, rng, white_ffts, initial_states, phi; store_paths)

Given an input trial C, perform the following Monte-Carlo map:

  D_eta(tau) = g^2 C(tau)
      -> draw eta paths
      -> solve the scalar nonlinear ODE
      -> C_output(tau) = E[phi(x(t)) phi(x(t+tau))].

Returns C_output, the ensemble mean of phi(x), and optionally a few paths.
"""
function evaluate_dmft_map(
    C::AbstractVector{<:Real},
    p::DMFTParams,
    rng::AbstractRNG,
    white_ffts::Union{Nothing, Matrix{ComplexF64}},
    initial_states::AbstractVector{<:Real},
    phi;
    store_paths::Bool = false,
)
    n = p.n_time
    length(C) == n || throw(DimensionMismatch(
        "C has length $(length(C)), but n_time=$n"
    ))
    length(initial_states) == p.n_paths || throw(DimensionMismatch(
        "initial_states has the wrong number of paths"
    ))

    # DMFT-2: D_eta = g^2 C.
    D_eta = (p.g^2) .* Float64.(C)

    # The FFT entries are the eigenvalues of the circulant covariance matrix.
    # They should be nonnegative for a valid covariance.  Clipping only protects
    # against tiny negative roundoff errors.
    covariance_eigenvalues = real.(fft(D_eta))
    covariance_eigenvalues .= max.(covariance_eigenvalues, 0.0)
    sqrt_eigenvalues = sqrt.(covariance_eigenvalues)

    C_output = zeros(Float64, n)
    mean_rate_sum = 0.0

    n_store = store_paths ? p.n_store_paths : 0
    sample_eta = Matrix{Float64}(undef, n, n_store)
    sample_x = Matrix{Float64}(undef, n, n_store)
    sample_rate = Matrix{Float64}(undef, n, n_store)

    for path in 1:p.n_paths
        # Either reuse the same standard Gaussian Fourier coefficients on every
        # fixed-point iteration or draw a fresh white-noise path.
        white_fft = if white_ffts === nothing
            fft(randn(rng, n))
        else
            @view white_ffts[:, path]
        end

        eta = colored_noise_from_white_fft(sqrt_eigenvalues, white_fft)
        x = integrate_effective_path(eta, p, initial_states[path], phi)
        rate = phi.(x)

        # Time-origin average for this path, then accumulated ensemble average.
        C_output .+= circular_autocorrelation(rate)
        mean_rate_sum += sum(rate)

        if path <= n_store
            sample_eta[:, path] .= eta
            sample_x[:, path] .= x
            sample_rate[:, path] .= rate
        end
    end

    C_output ./= p.n_paths
    C_output = project_to_circulant_covariance(C_output)
    mean_rate = mean_rate_sum / (p.n_paths * n)

    return C_output, mean_rate, sample_eta, sample_x, sample_rate
end


# ------------------------------------------------------------------------------
# Residuals and frequency grid
# ------------------------------------------------------------------------------

"""
    covariance_residual(C_output, C_input)

Compute RMS absolute and relative residuals over nonnegative lags 0 <= tau <= T/2.
The second half of a periodic covariance duplicates the corresponding negative
lags and is therefore omitted from the norm.
"""
function covariance_residual(
    C_output::AbstractVector{<:Real},
    C_input::AbstractVector{<:Real},
)
    n = length(C_input)
    length(C_output) == n || throw(DimensionMismatch("Covariance lengths differ"))

    n_nonnegative = div(n, 2) + 1
    output_half = @view C_output[1:n_nonnegative]
    input_half = @view C_input[1:n_nonnegative]

    absolute_rms = norm(output_half .- input_half) / sqrt(n_nonnegative)

    output_rms = norm(output_half) / sqrt(n_nonnegative)
    input_rms = norm(input_half) / sqrt(n_nonnegative)
    scale = max(output_rms, input_rms, 1.0e-14)

    relative_rms = absolute_rms / scale
    return relative_rms, absolute_rms
end


"""
    fft_angular_frequencies(n, dt)

Return angular frequencies in the ordering used by `fft` for even n.
"""
function fft_angular_frequencies(n::Int, dt::Float64)
    iseven(n) || throw(ArgumentError("n must be even"))

    integer_frequencies = vcat(
        collect(0:div(n, 2)),
        collect(-(div(n, 2) - 1):-1),
    )

    return (2.0 * pi / (n * dt)) .* Float64.(integer_frequencies)
end


# ------------------------------------------------------------------------------
# Main DMFT fixed-point solver
# ------------------------------------------------------------------------------

"""
    solve_dmft(p; initial_C=nothing, phi=tanh)

Solve the stationary dynamic mean-field equations by a damped Monte-Carlo
fixed-point iteration.

`initial_C` can be supplied to continue a solution from a nearby value of s or
g.  It must contain the full periodic covariance array with length p.n_time.
Continuation is particularly useful for following the dynamic solution into a
region where static fixed-point solutions coexist.

`phi` defaults to `tanh`.  A different scalar transfer function may be supplied,
but the same phi is then used in both the self-interaction and the random
network output, corresponding to

    dx_i/dt = -x_i + s*phi(x_i) + g*sum_j J_ij*phi(x_j).
"""
function solve_dmft(
    p::DMFTParams;
    initial_C::Union{Nothing, AbstractVector{<:Real}} = nothing,
    phi = tanh,
)
    validate_params(p)

    rng = MersenneTwister(p.seed)

    C = if initial_C === nothing
        initial_covariance(p)
    else
        length(initial_C) == p.n_time || throw(DimensionMismatch(
            "initial_C must have length n_time=$(p.n_time)"
        ))
        project_to_circulant_covariance(initial_C)
    end

    initial_states = make_initial_states(p, rng)

    white_ffts = if p.common_random_numbers
        make_white_noise_ffts(p, rng)
    else
        nothing
    end

    relative_history = Float64[]
    absolute_history = Float64[]
    step_history = Float64[]
    C0_history = Float64[]
    mean_rate_history = Float64[]

    converged = false
    iterations_completed = 0

    if p.verbose
        println("\nNonlinear self-interaction DMFT")
        println("--------------------------------")
        @printf("s = %.6g, g = %.6g\n", p.s, p.g)
        @printf("dt = %.6g, n_time = %d, T = %.6g\n",
                p.dt, p.n_time, p.dt * p.n_time)
        @printf("n_paths = %d, warmup_periods = %d, mixing = %.4g\n",
                p.n_paths, p.warmup_periods, p.mixing)
        println("Columns: iteration, C(0), relative residual, absolute residual, step RMS, mean tanh(x)")
    end

    for iteration in 1:p.max_iterations
        C_output, mean_rate, _, _, _ = evaluate_dmft_map(
            C,
            p,
            rng,
            white_ffts,
            initial_states,
            phi;
            store_paths = false,
        )

        relative_residual, absolute_residual = covariance_residual(C_output, C)

        # Damped Picard iteration of the functional fixed-point equation C=F[C].
        C_next = (1.0 - p.mixing) .* C .+ p.mixing .* C_output
        C_next = project_to_circulant_covariance(C_next)

        _, step_rms = covariance_residual(C_next, C)

        push!(relative_history, relative_residual)
        push!(absolute_history, absolute_residual)
        push!(step_history, step_rms)
        push!(C0_history, C_next[1])
        push!(mean_rate_history, mean_rate)

        iterations_completed = iteration
        C = C_next

        if p.verbose && (iteration == 1 || iteration % 5 == 0)
            @printf(
                "%6d  %11.6g  %12.5e  %12.5e  %12.5e  %+12.5e\n",
                iteration,
                C[1],
                relative_residual,
                absolute_residual,
                step_rms,
                mean_rate,
            )
        end

        if iteration >= p.min_iterations &&
           (relative_residual <= p.relative_tolerance ||
            absolute_residual <= p.absolute_tolerance)
            converged = true
            break
        end
    end

    # Evaluate the map once more at the returned covariance.  Besides providing
    # a final consistency check, this creates representative eta/x/rate traces.
    C_from_paths, final_mean_rate, sample_eta, sample_x, sample_rate =
        evaluate_dmft_map(
            C,
            p,
            rng,
            white_ffts,
            initial_states,
            phi;
            store_paths = true,
        )

    final_relative, final_absolute = covariance_residual(C_from_paths, C)

    n = p.n_time
    lags = p.dt .* Float64.(collect(0:div(n, 2)))
    angular_frequencies = fft_angular_frequencies(n, p.dt)

    # Multiplication by dt gives a Riemann-sum approximation to the continuous
    # Fourier transform of C(tau).  Tiny negative values are clipped.
    C_spectrum = p.dt .* max.(real.(fft(C)), 0.0)
    noise_spectrum = (p.g^2) .* C_spectrum
    noise_covariance = (p.g^2) .* C

    if p.verbose
        println()
        if converged
            println("Converged according to the requested tolerance.")
        else
            println("Reached max_iterations before satisfying the tolerance.")
        end
        @printf("Final C(0)               = %.8g\n", C[1])
        @printf("Final relative residual  = %.6e\n", final_relative)
        @printf("Final absolute residual  = %.6e\n", final_absolute)
        @printf("Final mean tanh(x)        = %+.6e\n", final_mean_rate)
        println("A nonzero mean much larger than Monte-Carlo error can indicate")
        println("asymmetric basin selection; use :bimodal or more paths for the")
        println("symmetric zero-mean solution.")
    end

    return DMFTResult(
        p,
        C,
        C_from_paths,
        noise_covariance,
        lags,
        angular_frequencies,
        C_spectrum,
        noise_spectrum,
        relative_history,
        absolute_history,
        step_history,
        C0_history,
        mean_rate_history,
        converged,
        iterations_completed,
        final_relative,
        final_absolute,
        final_mean_rate,
        sample_eta,
        sample_x,
        sample_rate,
    )
end


# ------------------------------------------------------------------------------
# Plain CSV output without an additional data-frame package
# ------------------------------------------------------------------------------

"""
    save_result(result; prefix="dmft")

Write four human-readable CSV files:

  * `<prefix>_correlation.csv`
  * `<prefix>_spectrum.csv`
  * `<prefix>_convergence.csv`
  * `<prefix>_sample_paths.csv`

The function returns a named tuple containing the filenames.
"""
function save_result(result::DMFTResult; prefix::AbstractString = "dmft")
    output_directory = dirname(prefix)
    output_directory != "." && mkpath(output_directory)

    correlation_file = prefix * "_correlation.csv"
    spectrum_file = prefix * "_spectrum.csv"
    convergence_file = prefix * "_convergence.csv"
    paths_file = prefix * "_sample_paths.csv"

    # Only nonnegative lags up to T/2 are written; the remaining half is the
    # even, periodic mirror image.
    n_nonnegative = length(result.lags)

    open(correlation_file, "w") do io
        println(io, "lag,C_iterate,C_from_paths,noise_covariance_g2C")
        for k in 1:n_nonnegative
            @printf(
                io,
                "%.16g,%.16g,%.16g,%.16g\n",
                result.lags[k],
                result.C[k],
                result.C_from_paths[k],
                result.noise_covariance[k],
            )
        end
    end

    open(spectrum_file, "w") do io
        println(io, "angular_frequency,C_spectrum,noise_spectrum")
        for k in eachindex(result.angular_frequencies)
            @printf(
                io,
                "%.16g,%.16g,%.16g\n",
                result.angular_frequencies[k],
                result.C_spectrum[k],
                result.noise_spectrum[k],
            )
        end
    end

    open(convergence_file, "w") do io
        println(io, "iteration,relative_residual,absolute_residual,step_rms,C0,mean_rate")
        for k in eachindex(result.relative_residual_history)
            @printf(
                io,
                "%d,%.16g,%.16g,%.16g,%.16g,%.16g\n",
                k,
                result.relative_residual_history[k],
                result.absolute_residual_history[k],
                result.step_history[k],
                result.C0_history[k],
                result.mean_rate_history[k],
            )
        end
    end

    open(paths_file, "w") do io
        n_store = size(result.sample_x, 2)
        header = String["time"]

        for path in 1:n_store
            push!(header, "eta_$path")
            push!(header, "x_$path")
            push!(header, "tanh_x_$path")
        end

        println(io, join(header, ","))

        for k in 1:result.params.n_time
            @printf(io, "%.16g", (k - 1) * result.params.dt)

            for path in 1:n_store
                @printf(
                    io,
                    ",%.16g,%.16g,%.16g",
                    result.sample_eta[k, path],
                    result.sample_x[k, path],
                    result.sample_rate[k, path],
                )
            end

            println(io)
        end
    end

    return (
        correlation = correlation_file,
        spectrum = spectrum_file,
        convergence = convergence_file,
        sample_paths = paths_file,
    )
end

end # module NonlinearSelfDMFT


# ==============================================================================
# Example run
# ==============================================================================
# This block is executed only when the file itself is run as a script.  It is
# skipped when the file is loaded with `include("nonlinear_self_dmft.jl")`.

if abspath(PROGRAM_FILE) == @__FILE__
    # Parameters close to those used in the paper's dynamic-DMFT plots.
    #
    # For a quick smoke test, try n_time=1024, n_paths=32, max_iterations=20.
    # For quantitative work, increase n_time and n_paths and verify convergence
    # with respect to dt, T=n_time*dt, the Monte-Carlo ensemble, and the seed.
    params = NonlinearSelfDMFT.DMFTParams(
        s = 1.2,
        g = 1.5,
        dt = 0.05,
        n_time = 4096,
        n_paths = 128,
        warmup_periods = 2,
        mixing = 0.20,
        max_iterations = 100,
        min_iterations = 10,
        relative_tolerance = 2.0e-3,
        absolute_tolerance = 5.0e-4,
        initial_x_mode = :auto,
        common_random_numbers = true,
        n_store_paths = 4,
        seed = 12345,
        verbose = true,
    )

    result = NonlinearSelfDMFT.solve_dmft(params)
    files = NonlinearSelfDMFT.save_result(result; prefix = "nonlinear_self_dmft")

    println("\nWrote:")
    println("  ", files.correlation)
    println("  ", files.spectrum)
    println("  ", files.convergence)
    println("  ", files.sample_paths)

    # Optional plotting, after installing Plots.jl:
    #
    #   using Plots
    #   nlag = length(result.lags)
    #   plot(result.lags, result.C[1:nlag],
    #        xlabel="lag tau", ylabel="C(tau)", label="DMFT")
    #
    # Continuation in s can help follow the dynamic branch when static and
    # dynamic solutions coexist:
    #
    #   p2 = NonlinearSelfDMFT.DMFTParams(s=1.3, g=1.5, dt=params.dt,
    #                   n_time=params.n_time, n_paths=params.n_paths)
    #   result2 = NonlinearSelfDMFT.solve_dmft(p2; initial_C=result.C)
end
