# This file was generated, do not modify it. # hide
using Random, Statistics, Plots

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
        return -√(2log(2))
    end
end


function simulate_curve(N::Int, βs::Vector{Float64}; nsamples::Int=30, seed::Int=1)
    # rng = MersenneTwister(seed)
    means = Float64[]
    stderrs = Float64[]
    for β in βs
        vals = [rem_free_energy_sample(N, β) for _ in 1:nsamples]
        push!(means, mean(vals))
        push!(stderrs, std(vals) / √nsamples)
        println("β=$(round(β,digits=3))  f_sim=$(round(vals[1],digits=4))  f_th=$(round(rem_free_energy_theory(β),digits=4))")
    end
    return means, stderrs
end


function plot_rem_free_energy(; N=18, βmin = 0.25, βmax = 2.25, nβ=25,
                              nsamples=30,
                              seed=1,
                              outfile="rem_free_energy.png")
    βc = √(2log(2))
    βs = collect(range(βmin, βmax, length=nβ))

    f_sim, f_err = simulate_curve(N, βs; nsamples=nsamples, seed=seed)

    p = plot(rem_free_energy_theory, xlim=(βmin, βmax);
        label="Theory (N → ∞)",
        linewidth=3,
        xlabel="β",
        ylabel="free energy f(β)",
        title="REM: theory vs simulation (N=$N, samples=$nsamples)"
    )

    scatter!(p, βs, f_sim;
        yerror=f_err,
        label="Simulation",
        markersize=3
    )

    vline!(p, [βc]; label="βc = √(2 log 2)", linestyle=:dash)

    savefig(p, outfile)
    #return p
end

plot_rem_free_energy(N=16, nsamples=30, outfile="./_assets/images/blog/rem_free_energy.png")