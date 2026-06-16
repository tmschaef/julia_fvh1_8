# Julia port of ppmlr.f90

"""
    flatten0!(sweeps::Sweeps)

Set all flattening coefficients to zero.
"""
function flatten0!(sweeps::Sweeps)
    nmin = sweeps.nmin
    nmax = sweeps.nmax
    for n in (nmin - 4):(nmax + 4)
        sweeps.flat[n] = 0.0
    end
    return nothing
end

"""
    riemann0!(lmin, lmax, gamma, prgh, urgh, vrgh, plft, ulft, vlft, pmid, umid)

Simple Riemann average used by the Julia port.
"""
function riemann0!(lmin::Integer, lmax::Integer, gamma::Real,
                   prgh::AbstractVector{<:Real}, urgh::AbstractVector{<:Real},
                   vrgh::AbstractVector{<:Real}, plft::AbstractVector{<:Real},
                   ulft::AbstractVector{<:Real}, vlft::AbstractVector{<:Real},
                   pmid::AbstractVector{<:Real}, umid::AbstractVector{<:Real})
    smallp_local = 1.0e-25
    for l in lmin:lmax
        umid[l] = 0.5 * (urgh[l] + ulft[l])
        pmid[l] = 0.5 * (prgh[l] + plft[l])
        pmid[l] = max(smallp_local, pmid[l])
    end
    return nothing
end

"""
    volume!(sweeps::Sweeps)

Compute zone volumes and copy the current geometry into the remap arrays.
"""
function volume!(sweeps::Sweeps)
    nmin = sweeps.nmin
    nmax = sweeps.nmax
    sweeps.radius = 1.0
    for n in (nmin - 3):(nmax + 4)
        sweeps.dvol[n] = sweeps.dx[n]
        sweeps.dvol0[n] = sweeps.dx0[n]
    end
    return nothing
end

"""
    ppmlr!(sweeps, g, w)

Run the PPM/Lagrangian-remap workflow.
"""
function ppmlr!(sweeps::Sweeps, g::globalvars, w::PpmWork)
    nmin = sweeps.nmin
    nmax = sweeps.nmax
    para = w.para

    boundary!(sweeps)
    paraset!(para)
    volume!(sweeps)
    flatten0!(sweeps)

    dr = w.dr
    du = w.du
    dp = w.dp
    r6 = w.r6
    u6 = w.u6
    p6 = w.p6
    rl = w.rl
    ul = w.ul
    pl = w.pl
    rrgh = w.rrgh
    urgh = w.urgh
    prgh = w.prgh
    rlft = w.rlft
    ulft = w.ulft
    plft = w.plft
    umid = w.umid
    pmid = w.pmid

    parabola!(nmin - 4, nmax + 4, para, sweeps.p, dp, p6, pl, sweeps.flat, w)
    parabola!(nmin - 4, nmax + 4, para, sweeps.r, dr, r6, rl, sweeps.flat, w)
    parabola!(nmin - 4, nmax + 4, para, sweeps.u, du, u6, ul, sweeps.flat, w)

    states!(sweeps, g, pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh, w)
    riemann0!(nmin - 3, nmax + 4, gam, prgh, urgh, rrgh, plft, ulft, rlft, pmid, umid)

    evolve!(sweeps, g, umid, pmid, w)
    remap!(sweeps, w)

    return nothing
end
