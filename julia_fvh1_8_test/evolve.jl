# Julia port of evolve.f90
# Depends on the already-porting modules in vh1mods.jl.

include("vh1mods.jl")

using .Global
using .Sweeps
using .Sweepsize

"""
    evolve!(umid, pmid)

Lagrangian update for density, velocity, and total energy.
"""
function evolve!(umid::AbstractVector{<:Real}, pmid::AbstractVector{<:Real})
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    # Local work arrays
    amid = zeros(Float64, Sweepsize.maxsweep)
    uold = zeros(Float64, Sweepsize.maxsweep)
    xa1  = zeros(Float64, Sweepsize.maxsweep)
    dvol1 = zeros(Float64, Sweepsize.maxsweep)
    upmid = zeros(Float64, Sweepsize.maxsweep)
    dm    = zeros(Float64, Sweepsize.maxsweep)
    dtbdm = zeros(Float64, Sweepsize.maxsweep)
    xa2   = zeros(Float64, Sweepsize.maxsweep)
    xa3   = zeros(Float64, Sweepsize.maxsweep)

    # Grid position evolution
    for n in (nmin - 3):(nmax + 4)
        dm[n]    = Sweeps.r[n] * Sweeps.dvol[n]
        dtbdm[n] = Global.dt / dm[n]
        xa1[n]   = Sweeps.xa[n]
        dvol1[n] = Sweeps.dvol[n]
        Sweeps.xa[n]    = Sweeps.xa[n] + Global.dt * umid[n] / Sweeps.radius
        upmid[n] = umid[n] * pmid[n]
    end

    xa1[nmin - 4] = Sweeps.xa[nmin - 4]
    xa1[nmax + 5] = Sweeps.xa[nmax + 5]

    for n in (nmin - 4):(nmax + 5)
        xa2[n] = xa1[n] + 0.5 * Sweeps.dx[n]
        Sweeps.dx[n]  = Sweeps.xa[n + 1] - Sweeps.xa[n]
        xa3[n] = Sweeps.xa[n] + 0.5 * Sweeps.dx[n]
    end

    # Calculate dvolume and average area
    for n in (nmin - 3):(nmax + 4)
        Sweeps.dvol[n] = Sweeps.dx[n]
        amid[n] = 1.0
    end

    for n in (nmin - 3):(nmax + 3)
        # Density evolution. Lagrangian code, only change volume.
        Sweeps.r[n] = Sweeps.r[n] * (dvol1[n] / Sweeps.dvol[n])
        Sweeps.r[n] = max(Sweeps.r[n], Global.smallr)

        # Velocity evolution due to pressure acceleration and forces.
        uold[n] = Sweeps.u[n]
        Sweeps.u[n] = Sweeps.u[n] - dtbdm[n] * (pmid[n + 1] - pmid[n]) * 0.5 * (amid[n + 1] + amid[n])

        # Total energy evolution
        Sweeps.e[n] = Sweeps.e[n] - dtbdm[n] * (amid[n + 1] * upmid[n + 1] - amid[n] * upmid[n])
        Sweeps.q[n] = Sweeps.e[n] - 0.5 * (Sweeps.u[n]^2 + Sweeps.v[n]^2 + Sweeps.w[n]^2)
        Sweeps.p[n] = max(Sweeps.r[n] * Sweeps.q[n] * Global.gamm, Global.smallp)
    end

    return nothing
end
