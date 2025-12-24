# This file was generated, do not modify it. # hide
using Random, Statistics, Plots, LaTeXStrings

"""
    Robust implementation of logsumexp
"""
logsumexp(x) = (m = maximum(x); m + log(sum(exp.(x .- m))))

"""
    One-sample REM free energy
"""
function rem_free_energy_sample(N::Int, β::Float64)
    M = 1 << N
    σ = sqrt(N)
    E = σ .* randn(M)
    logZ = logsumexp(- β .* E)
    return - logZ / (β * N)
end


"""
    Averaged REM free energy in the limit N → ∞
"""
function rem_free_energy_theory(β::Float64)
    βc = √(2log(2))
    if β <= βc
        return - log(2) / β - β / 2
    else
        return - √(2log(2))
    end
end

"""
    Computes the average over `nsamples` of the REM free energy
    for the values of β contained in the vector `βs`
"""
function simulate_curve(N::Int, βs::Vector{Float64}; nsamples::Int=30, seed::Int=23)
    if seed > 0
        Random.seed!(seed);
    end
    means = Float64[]
    stderrs = Float64[]
    for β in βs
        vals = [rem_free_energy_sample(N, β) for _ in 1:nsamples]
        push!(means, mean(vals))
        push!(stderrs, std(vals) / √nsamples)
    end
    return means, stderrs
end

"""
    Generates a plot comparing theory vs simulations
"""
function plot_rem_free_energy(; N=18, βmin = 0.25, βmax = 2.25, nβ=25,
                              nsamples=30,
                              seed=-23,
                              outfile="rem_free_energy.png")
    βc = √(2log(2))
    βs = collect(range(βmin, βmax, length=nβ))

    f_sim, f_err = simulate_curve(N, βs; nsamples=nsamples, seed=seed)

    plot(rem_free_energy_theory, xlim=(βmin, βmax);
        label="Theory (N → ∞)",
        linewidth=3, color=:red,
        xlabel="β", ylabel="free energy f(β)",
    )

    scatter!(βs, f_sim;
        yerror=f_err, color=:blue, markerstrokecolor=:blue,
        label="Simulation",
        markersize=3.5
    )

    vline!([βc]; label=L"β_c"*"= √(2 log 2)", linestyle=:dash, linewidth=3, color=:green)

    savefig(outfile)
end

plot_rem_free_energy(N=16, nsamples=30, outfile="./_assets/images/blog/rem_free_energy.png")