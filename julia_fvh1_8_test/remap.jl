# Julia port of remap.f90
# Remaps conserved quantities from the Lagrangian grid to the Eulerian grid.

include("vh1mods.jl")
include("parabola.jl")

using .Global
using .Sweeps
using .Sweepsize

"""
    remap!()

Apply the remap update from the Lagrangian grid to the Eulerian grid.
"""
function remap!()
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    third = 1.0 / 3.0
    twothd = 2.0 / 3.0
    fourthd = 4.0 / 3.0

    para::Vector{Float64} = zeros(Float64, 5)
    
    du = zeros(Float64, Sweepsize.maxsweep)
    ul = zeros(Float64, Sweepsize.maxsweep)
    u6 = zeros(Float64, Sweepsize.maxsweep)
    dv = zeros(Float64, Sweepsize.maxsweep)
    vl = zeros(Float64, Sweepsize.maxsweep)
    v6 = zeros(Float64, Sweepsize.maxsweep)
    dw = zeros(Float64, Sweepsize.maxsweep)
    wl = zeros(Float64, Sweepsize.maxsweep)
    w6 = zeros(Float64, Sweepsize.maxsweep)
    de = zeros(Float64, Sweepsize.maxsweep)
    el = zeros(Float64, Sweepsize.maxsweep)
    e6 = zeros(Float64, Sweepsize.maxsweep)

    dq = zeros(Float64, Sweepsize.maxsweep)
    ql = zeros(Float64, Sweepsize.maxsweep)
    q6 = zeros(Float64, Sweepsize.maxsweep)
    dr = zeros(Float64, Sweepsize.maxsweep)
    rl = zeros(Float64, Sweepsize.maxsweep)
    r6 = zeros(Float64, Sweepsize.maxsweep)
    dm = zeros(Float64, Sweepsize.maxsweep)
    dm0 = zeros(Float64, Sweepsize.maxsweep)
    delta = zeros(Float64, Sweepsize.maxsweep)

    fluxr = zeros(Float64, Sweepsize.maxsweep)
    fluxu = zeros(Float64, Sweepsize.maxsweep)
    fluxv = zeros(Float64, Sweepsize.maxsweep)
    fluxw = zeros(Float64, Sweepsize.maxsweep)
    fluxe = zeros(Float64, Sweepsize.maxsweep)
    fluxq = zeros(Float64, Sweepsize.maxsweep)

    paraset!(para)
    parabola!(nmin - 1, nmax + 1, para, Sweeps.r, dr, r6, rl, Sweeps.flat)
    parabola!(nmin - 1, nmax + 1, para, Sweeps.u, du, u6, ul, Sweeps.flat)
    parabola!(nmin - 1, nmax + 1, para, Sweeps.v, dv, v6, vl, Sweeps.flat)
    parabola!(nmin - 1, nmax + 1, para, Sweeps.w, dw, w6, wl, Sweeps.flat)
    parabola!(nmin - 1, nmax + 1, para, Sweeps.q, dq, q6, ql, Sweeps.flat)
    parabola!(nmin - 1, nmax + 1, para, Sweeps.e, de, e6, el, Sweeps.flat)

    for n in nmin:(nmax + 1)
        delta[n] = Sweeps.xa[n] - Sweeps.xa0[n]
    end

    for n in nmin:(nmax + 1)
        deltx = Sweeps.xa[n] - Sweeps.xa0[n]
        if deltx >= 0.0
            nn = n - 1
            fractn = 0.5 * deltx / Sweeps.dx[nn]
            fractn2 = 1.0 - fourthd * fractn
            fluxr[n] = (rl[nn] + dr[nn] - fractn * (dr[nn] - fractn2 * r6[nn])) * delta[n]
            fluxu[n] = (ul[nn] + du[nn] - fractn * (du[nn] - fractn2 * u6[nn])) * fluxr[n]
            fluxv[n] = (vl[nn] + dv[nn] - fractn * (dv[nn] - fractn2 * v6[nn])) * fluxr[n]
            fluxw[n] = (wl[nn] + dw[nn] - fractn * (dw[nn] - fractn2 * w6[nn])) * fluxr[n]
            fluxe[n] = (el[nn] + de[nn] - fractn * (de[nn] - fractn2 * e6[nn])) * fluxr[n]
            fluxq[n] = (ql[nn] + dq[nn] - fractn * (dq[nn] - fractn2 * q6[nn])) * fluxr[n]
        else
            fractn = 0.5 * deltx / Sweeps.dx[n]
            fractn2 = 1.0 + fourthd * fractn
            fluxr[n] = (rl[n] - fractn * (dr[n] + fractn2 * r6[n])) * delta[n]
            fluxu[n] = (ul[n] - fractn * (du[n] + fractn2 * u6[n])) * fluxr[n]
            fluxv[n] = (vl[n] - fractn * (dv[n] + fractn2 * v6[n])) * fluxr[n]
            fluxw[n] = (wl[n] - fractn * (dw[n] + fractn2 * w6[n])) * fluxr[n]
            fluxe[n] = (el[n] - fractn * (de[n] + fractn2 * e6[n])) * fluxr[n]
            fluxq[n] = (ql[n] - fractn * (dq[n] + fractn2 * q6[n])) * fluxr[n]
        end
    end

    for n in nmin:nmax
        dm[n] = Sweeps.r[n] * Sweeps.dvol[n]
        dm0[n] = dm[n] + fluxr[n] - fluxr[n + 1]
        Sweeps.r[n] = dm0[n] / Sweeps.dvol0[n]
        Sweeps.r[n] = max(Global.smallr, Sweeps.r[n])
        dm0[n] = 1.0 / (Sweeps.r[n] * Sweeps.dvol0[n])
        Sweeps.u[n] = (Sweeps.u[n] * dm[n] + fluxu[n] - fluxu[n + 1]) * dm0[n]
        Sweeps.v[n] = (Sweeps.v[n] * dm[n] + fluxv[n] - fluxv[n + 1]) * dm0[n]
        Sweeps.w[n] = (Sweeps.w[n] * dm[n] + fluxw[n] - fluxw[n + 1]) * dm0[n]
        Sweeps.e[n] = (Sweeps.e[n] * dm[n] + fluxe[n] - fluxe[n + 1]) * dm0[n]
        Sweeps.q[n] = (Sweeps.q[n] * dm[n] + fluxq[n] - fluxq[n + 1]) * dm0[n]

        ekin = 0.5 * (Sweeps.u[n]^2 + Sweeps.v[n]^2 + Sweeps.w[n]^2)
        if ekin / Sweeps.q[n] < 100.0
            Sweeps.q[n] = Sweeps.e[n] - ekin
        end

        Sweeps.p[n] = Global.gamm * Sweeps.r[n] * Sweeps.q[n]
        Sweeps.p[n] = max(Global.smallp, Sweeps.p[n])
    end

    return nothing
end
