# This file was generated, do not modify it. # hide
using Random, Statistics
using Plots

# ---------- Numerics ----------
logsumexp(x) = (m = maximum(x); m + log(sum(exp.(x .- m))))

"""
One-sample REM free energy per spin:
Eα ~ Normal(0, N), α=1..2^N
fN(β) = -(1/(βN)) log Σ exp(-β Eα)
"""
function rem_free_energy_per_spin(N::Int, β::Float64; rng=Random.default_rng())
    M = 1 << N
    σ = sqrt(N)
    E = σ .* randn(rng, M)
    logZ = logsumexp(-β .* E)
    return -(1 / (β * N)) * logZ
end

# ---------- Theory (N→∞) ----------
function rem_free_energy_theory(β::Float64)
    βc = sqrt(2 * log(2))
    if β <= βc
        return -(log(2)/β) - (β/2)
    else
        return -sqrt(2 * log(2))
    end
end

# ---------- Experiment: average over disorder ----------
function simulate_curve(N::Int, betas::Vector{Float64}; nsamples::Int=30, seed::Int=1)
    rng = MersenneTwister(seed)
    means = Float64[]
    stderrs = Float64[]
    for β in betas
        vals = [rem_free_energy_per_spin(N, β; rng=rng) for _ in 1:nsamples]
        push!(means, mean(vals))
        push!(stderrs, std(vals) / sqrt(nsamples))
    end
    return means, stderrs
end

# ---------- Plot ----------
function plot_rem_free_energy(; N=18,
                              betas=collect(range(0.25, 2.25, length=25)),
                              nsamples=30,
                              seed=1,
                              outfile="rem_free_energy.png")
    βc = sqrt(2 * log(2))
    f_th = [rem_free_energy_theory(β) for β in betas]
    f_sim, f_err = simulate_curve(N, betas; nsamples=nsamples, seed=seed)

    p = plot(betas, f_th;
        label="Theory (N→∞)",
        linewidth=3,
        xlabel="β",
        ylabel="free energy f(β)",
        title="REM: theory vs simulation (N=$N, samples=$nsamples)"
    )

    scatter!(p, betas, f_sim;
        yerror=f_err,
        label="Simulation (mean ± s.e.)",
        markersize=4
    )

    vline!(p, [βc]; label="βc = √(2 log 2)", linestyle=:dash)

    savefig(p, outfile)
    return p
end

plot_rem_free_energy(N=16, nsamples=30, outfile="./_assets/images/blog/rem_free_energy.png")