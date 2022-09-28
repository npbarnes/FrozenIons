module FrozenIons

using Unitful
using StatsBase
using CSV
using DifferentialEquations
using PyPlot: plt
using QuadGK

export plot_overlay, simulate

const _rocketvelocity = 2.5u"km/s"
const _N0 = 2.631e24
const _τp = 28u"s"
const _τb = 5u"s"
const _FWHMrate = 3u"km/s"
const _dσ = _FWHMrate/2√(2log(2))
const _vth = 1.8u"km/s"

P_thermal(v; vth=_vth) = 4π * (1/(π*vth^2))^(3/2) * v^2 * exp(-v^2/vth^2)
P_CRRES(v) = error("Not Implemented")

"""Nn is zero until t=0 when N0 neutrals become available. The value is constant for t>0."""
Nn_constant(t; N0=_N0) = t <= 0u"s" ? zero(N0) : N0
"""At t=0 N0 neutrals are made available, but photoionize with time constant τp"""
Nn_photo(t; N0=_N0, τp=_τp, kwargs...) =  t <= 0u"s" ? zero(N0) : N0*exp(-t/τp)
"""Neutrals are made available gradually with time constant τb, but is simultaineously depleted by photoionization with time constant τp"""
Nn_slowburn(t; N0=_N0, τp=_τp, τb=_τb) = t <= 0u"s" ? zero(N0) : N0/(1-τb/τp) * (exp(-t/τp) - exp(-t/τb))

"""Density of an instant localized release."""
n_instant(t,r; P=P_thermal, Nn=Nn_photo) = Nn(t)/(4π*r^2*t) * P(r/t)

"""Gaussian function but it's zero when σ is zero. The true limit would be a delta function."""
gauss(x; μ=zero(x), σ=oneunit(x)) = σ≈zero(σ) ? zero(inv(σ)) : 1/(σ*√(2π)) * exp(-1/2 * ((x-μ)/σ)^2)

simulate_slowburn(tspan, xs; Nn, kwargs...) = simulate(tspan, xs; Nn=Nn_slowburn, kwargs...)
function simulate(tspan, xs; Nn, N0=_N0, τp=_τp, τb=_τb, v=_rocketvelocity, dσ=_dσ, shape=gauss, toffset=0u"s")
    p = (;xs, N0, τp, τb, v, dσ, shape)
    function f(u,p,t)
        tt = t - toffset
        Nn(tt; p.N0, p.τp, p.τb)/p.τp * p.shape.(p.xs; μ=p.v*t, σ=p.dσ*tt)
    end
    prob = ODEProblem(f, zeros(typeof(inv(first(p.xs))), length(p.xs)), tspan, p)
    solve(prob, Tsit5(), tstops=toffset), p
end

function simulate_2d(tspan)
    qx = range(-25u"km", 250u"km", length=200)
    qy = range(-20u"km", 20u"km", length=100)
    u0 = zeros(typeof(inv(first(qx)*first(qy))), length(qx), length(qy))

    ñ(t,x,y; Nn, dσ) = Nn(t) * gauss(x; μ=0u"km", σ=dσ*t) * gauss(y; μ=0u"km", σ=dσ*t)

    RHS(u,p,t) = [ñ(t, x .- p.vx*t, y .- p.vy*t; p.Nn, p.dσ)/p.τp for x in qx, y in qy]

    p = (;qx, qy, vx=2.5u"km/s", vy=0u"km/s", τp=_τp, dσ=_dσ, Nn=Nn_slowburn)
    prob = ODEProblem(RHS, u0, tspan, p)
    solve(prob, Tsit5()), p
end

density_function(Nn, P) =  Nn(t)/(4π*r^2*t) * P(r/t)
function simulate_3d_mc(tspan, n_n)
    cloud_velocity = [2.5, 0, 0]u"km/s"
    cloud_center(t) = cloud_velocity .* t

    ts = range(0u"s", 100u"s", step=0.1u"s")
    dt = step(ts)

    for t in ts
        new_particles = [dσ*t*randn(3) + cloud_center(t) for _ in 1:round(Int, α*Nn(t))]
    end

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

plot_overlay(; Nn) = plot_overlay(plt.subplots()...; Nn)
function plot_overlay(fig, ax; Nn)#; v, τp, τb, R1overR2=1, R2_delay=50u"s")
    plot_donsions_normed(fig, ax)

    tspan = (0.0u"s", 1000.0u"s")
    xs = range(-50u"km", 500u"km", length=10000)
    # Solve release 1
    sol1, p1 = simulate(tspan, xs; Nn, v=10u"km/s")
    # Solve release 2
    sol2, p2 = simulate(tspan, xs; Nn, v=10u"km/s", toffset=50u"s")
    # add solutions together
    ssol = sol1[end] + sol2[end]
    # convert distance to time
    # plot
    ax.plot(ustrip.(u"s", p1.xs ./ p1.v), ssol./maximum(ssol))

    ax.set_xlim([-100,200])
    ax.set_xlabel("Time along trajectory (s)")
    ax.set_ylabel("Arbitrary (normalized) quantity")

    fig, ax
end

end # module
