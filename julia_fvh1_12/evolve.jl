# Julia port of evolve.f90

"""
    evolve!(sweeps, g, umid, pmid, w)

Lagrangian update for density, velocity, and total energy.
"""
function evolve!(sweeps::Sweeps, g::globalvars,
                 umid::AbstractVector{<:Real}, pmid::AbstractVector{<:Real},
                 w::PpmWork)
    nmin = sweeps.nmin
    nmax = sweeps.nmax

    amid = w.amid
    uold = w.uold
    xa1 = w.xa1
    dvol1 = w.dvol1
    upmid = w.upmid
    dm = w.dm
    dtbdm = w.dtbdm
    xa2 = w.xa2
    xa3 = w.xa3

    for n in (nmin - 3):(nmax + 4)
        dm[n]    = sweeps.r[n] * sweeps.dvol[n]
        dtbdm[n] = g.dt / dm[n]
        xa1[n]   = sweeps.xa[n]
        dvol1[n] = sweeps.dvol[n]
        sweeps.xa[n]    = sweeps.xa[n] + g.dt * umid[n] / sweeps.radius
        upmid[n] = umid[n] * pmid[n]
    end

    xa1[nmin - 4] = sweeps.xa[nmin - 4]
    xa1[nmax + 5] = sweeps.xa[nmax + 5]

    for n in (nmin - 4):(nmax + 5)
        xa2[n] = xa1[n] + 0.5 * sweeps.dx[n]
        sweeps.dx[n]  = sweeps.xa[n + 1] - sweeps.xa[n]
        xa3[n] = sweeps.xa[n] + 0.5 * sweeps.dx[n]
    end

    for n in (nmin - 3):(nmax + 4)
        sweeps.dvol[n] = sweeps.dx[n]
        amid[n] = 1.0
    end

    for n in (nmin - 3):(nmax + 3)
        sweeps.r[n] = sweeps.r[n] * (dvol1[n] / sweeps.dvol[n])
        sweeps.r[n] = max(sweeps.r[n], smallr)

        uold[n] = sweeps.u[n]
        sweeps.u[n] = sweeps.u[n] - dtbdm[n] * (pmid[n + 1] - pmid[n]) * 0.5 * (amid[n + 1] + amid[n])

        sweeps.e[n] = sweeps.e[n] - dtbdm[n] * (amid[n + 1] * upmid[n + 1] - amid[n] * upmid[n])
        sweeps.q[n] = sweeps.e[n] - 0.5 * (sweeps.u[n]^2 + sweeps.v[n]^2 + sweeps.w[n]^2)
        sweeps.p[n] = max(sweeps.r[n] * sweeps.q[n] * gamm, smallp)
    end

    return nothing
end
