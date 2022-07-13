module FrozenIons

using Unitful
using StatsBase
using CSV
using DifferentialEquations
using PyPlot: plt

export plot_overlay, Nn_constant, Nn_instant, Nn_slowburn, ionizationratedensity_1d

const _rocketvelocity = 2.5u"km/s"
const _samplescale = 10000/1e22
const _N0 = 2.631e24
const _τp = 28u"s"
const _τb = 5u"s"
const _FWHMrate = 3u"km/s"

σ(t; FWHMrate=_FWHMrate) = FWHMrate/2√(2log(2)) * t
dNidt(t; N0, τp, τb) = N0*(exp(-t/τp) - exp(-t/τb))/(τp - τb)

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

function sample1d(ts; v, τp, τb, burst_delay=0.0u"s", FWHMrate=_FWHMrate, N0=_N0, scale=_samplescale)
    dt = step(ts)
    result = []
    for t in ts
        n = floor(Int, scale*dt*dNidt(t; N0, τp, τb))
        xs = σ(t; FWHMrate).*randn(n) .+ v*(t + burst_delay)
        append!(result, xs)
    end
    result
end

function bincenters(edges)
    @. (edges[2:end] + edges[1:end-1]) / 2
end

plot_donsions_normed() = plot_donsions_normed(plt.subplots()...)
function plot_donsions_normed(fig, ax)
    f = donsions()
    T = f.time .* u"s"
    I_normed = f.value ./ maximum(f.value)

    ax.plot(ustrip.(u"s", T), I_normed)
end

plot_overlay(; kwargs...) = plot_overlay(plt.subplots()...; kwargs...)
function plot_overlay(fig, ax; v, τp, τb, R1overR2=1, R2_delay=50u"s")
    plot_donsions_normed(fig, ax)

    ts = range(0u"s", 1000u"s", length=4000+1)[2:end] # don't include zero
    samples = sample1d(ts; v, τp, τb, N0=R1overR2*_N0, scale=5*_samplescale) # Release 1
    append!(samples, sample1d(ts; v, τp, τb, scale=5*_samplescale, burst_delay=R2_delay)) # Release 2
    h = fit(Histogram, ustrip.(u"km", samples); nbins=8000)
    h_normed = h.weights ./ maximum(h.weights)
    bc = bincenters(only(h.edges))
    bc_astime = bc ./ ustrip(u"km/s", v)

    ax.plot(bc_astime, h_normed)

    ax.set_xlim([-100,200])
    ax.set_xlabel("Time along trajectory (s)")
    ax.set_ylabel("Arbitrary (normalized) quantity")

    fig, ax
end

Nn_constant(t; N0=_N0) = t <= 0u"s" ? zero(N0) : N0
Nn_instant(t; N0=_N0, τp=_τp) =  t <= 0u"s" ? zero(N0) : N0*exp(-t/τp)
Nn_slowburn(t; N0=_N0, τp=_τp, τb=_τb) = t <= 0u"s" ? zero(N0) : N0/(1-τb/τp) * (exp(-t/τp) - exp(-t/τb))

end # module
