# Julia port of ppmlr.f90
# Executes the 1D PPM/Lagrangian remap workflow using the translated helpers.

include("vh1mods.jl")
include("parabola.jl")
include("boundary.jl")
include("states.jl")
include("evolve.jl")
include("remap.jl")

using .Global
using .Sweeps
using .Sweepsize

"""
    flatten0!()

Set all flattening coefficients to zero.
"""
function flatten0!()
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    for n in (nmin - 4):(nmax + 4)
        Sweeps.flat[n] = 0.0
    end
    return nothing
end

"""
    riemann0!(lmin::Integer, lmax::Integer, gamma::Real,
              prgh, urgh, vrgh, plft, ulft, vlft, pmid, umid)

Simple Riemann average used by the Julia port.
"""
function riemann0!(lmin::Integer, lmax::Integer, gamma::Real,
                   prgh::AbstractVector{<:Real}, urgh::AbstractVector{<:Real},
                   vrgh::AbstractVector{<:Real}, plft::AbstractVector{<:Real},
                   ulft::AbstractVector{<:Real}, vlft::AbstractVector{<:Real},
                   pmid::AbstractVector{<:Real}, umid::AbstractVector{<:Real})
    smallp = 1.0e-25
    for l in lmin:lmax
        umid[l] = 0.5 * (urgh[l] + ulft[l])
        pmid[l] = 0.5 * (prgh[l] + plft[l])
        pmid[l] = max(smallp, pmid[l])
    end
    return nothing
end

"""
    volume!()

Compute zone volumes and copy the current geometry into the remap arrays.
"""
function volume!()
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    Sweeps.radius = 1.0
    for n in (nmin - 3):(nmax + 4)
        Sweeps.dvol[n] = Sweeps.dx[n]
        Sweeps.dvol0[n] = Sweeps.dx0[n]
    end
    return nothing
end

"""
    ppmlr!()

Run the PPM/Lagrangian-remap workflow.
"""
function ppmlr!()
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    para::Vector{Float64} = zeros(Float64, 5)
    
    boundary!()
    paraset!(para)
    volume!()
    flatten0!()

    dr = zeros(Float64, Sweepsize.maxsweep)
    du = zeros(Float64, Sweepsize.maxsweep)
    dp = zeros(Float64, Sweepsize.maxsweep)
    r6 = zeros(Float64, Sweepsize.maxsweep)
    u6 = zeros(Float64, Sweepsize.maxsweep)
    p6 = zeros(Float64, Sweepsize.maxsweep)
    rl = zeros(Float64, Sweepsize.maxsweep)
    ul = zeros(Float64, Sweepsize.maxsweep)
    pl = zeros(Float64, Sweepsize.maxsweep)

    rrgh = zeros(Float64, Sweepsize.maxsweep)
    urgh = zeros(Float64, Sweepsize.maxsweep)
    prgh = zeros(Float64, Sweepsize.maxsweep)
    rlft = zeros(Float64, Sweepsize.maxsweep)
    ulft = zeros(Float64, Sweepsize.maxsweep)
    plft = zeros(Float64, Sweepsize.maxsweep)
    umid = zeros(Float64, Sweepsize.maxsweep)
    pmid = zeros(Float64, Sweepsize.maxsweep)

    parabola!(nmin - 4, nmax + 4, para, Sweeps.p, dp, p6, pl, Sweeps.flat)
    parabola!(nmin - 4, nmax + 4, para, Sweeps.r, dr, r6, rl, Sweeps.flat)
    parabola!(nmin - 4, nmax + 4, para, Sweeps.u, du, u6, ul, Sweeps.flat)

    states!(pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh)
    riemann0!(nmin - 3, nmax + 4, Global.gam, prgh, urgh, rrgh, plft, ulft, rlft, pmid, umid)

    evolve!(umid, pmid)
    remap!()

    return nothing
end
