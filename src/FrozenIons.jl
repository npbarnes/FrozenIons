module FrozenIons

using Unitful
using StatsBase
using CSV
using DifferentialEquations
using PyPlot: plt

export plot_overlay

const _rocketvelocity = 2.5u"km/s"
const _N0 = 2.631e24
const _τp = 28u"s"
const _τb = 5u"s"
const _FWHMrate = 3u"km/s"
const _dσ = _FWHMrate/2√(2log(2))

Nn_constant(t; N0=_N0) = t <= 0u"s" ? zero(N0) : N0
Nn_instant(t; N0=_N0, τp=_τp) =  t <= 0u"s" ? zero(N0) : N0*exp(-t/τp)
Nn_slowburn(t; N0=_N0, τp=_τp, τb=_τb) = t <= 0u"s" ? zero(N0) : N0/(1-τb/τp) * (exp(-t/τp) - exp(-t/τb))

"""Gaussian function but it's zero when σ is zero. The true limit would be a delta function."""
gauss(x; μ=0, σ=1) = σ≈zero(σ) ? zero(inv(σ)) : 1/(σ*√(2π)) * exp(-1/2 * ((x-μ)/σ)^2)

function simulate_slowburn(tspan, xs; N0=_N0, τp=_τp, τb=_τb, v=_rocketvelocity, dσ=_dσ, shape=gauss)
    p = (;xs, N0, τp, τb, v, dσ, shape)
    f(u,p,t) = Nn_slowburn(t; p.N0, p.τp, p.τb)/p.τp * p.shape.(p.xs; μ=p.v*t, σ=p.dσ*t)
    prob = ODEProblem(f, zeros(typeof(inv(first(p.xs))), length(p.xs)), tspan, p)
    solve(prob, Tsit5()), p
end

function donsions()
    CSV.File("/home/nathan/Dropbox/AGU posters/CEDAR 2022/KX_BER_Z6_50mm_Ba_ion_profile.txt";
        header=["time","value"],
        delim=",",
        ignorerepeated=true,
        types=Float64
    )
end

function donsions_R1()
    f = donsions()
    (time=f.time[f.time .<= 45], value=f.value[f.time .<= 45])
end

plot_donsions_normed() = plot_donsions_normed(plt.subplots()...)
function plot_donsions_normed(fig, ax)
    f = donsions()
    ax.plot(f.time, f.value ./ maximum(f.value))
end

plot_overlay(; kwargs...) = plot_overlay(plt.subplots()...; kwargs...)
function plot_overlay(fig, ax)#; v, τp, τb, R1overR2=1, R2_delay=50u"s")
    plot_donsions_normed(fig, ax)

    tspan = (0.0u"s", 1000.0u"s")
    xs = range(-50u"km", 500u"km", length=10000)
    # Solve release 1
    sol, p = simulate_slowburn(tspan, xs)
    # Solve release 2
    # add solutions together
    # convert distance to time
    # plot
    ax.plot(ustrip.(u"s", p.xs ./ p.v), sol[end]./maximum(sol[end]))

    ax.set_xlim([-100,200])
    ax.set_xlabel("Time along trajectory (s)")
    ax.set_ylabel("Arbitrary (normalized) quantity")

    fig, ax
end

end # module
