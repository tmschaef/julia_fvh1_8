# Julia port of remap.f90

"""
    remap!(sweeps, w)

Apply the remap update from the Lagrangian grid to the Eulerian grid.
"""
function remap!(sweeps::Sweeps, w::PpmWork)
    nmin = sweeps.nmin
    nmax = sweeps.nmax
    fourthd = 4.0 / 3.0

    para = w.para
    du = w.du
    ul = w.ul
    u6 = w.u6
    dv = w.dv
    vl = w.vl
    v6 = w.v6
    dw = w.dw
    wl = w.wl
    w6 = w.w6
    de = w.de
    el = w.el
    e6 = w.e6
    dq = w.dq
    ql = w.ql
    q6 = w.q6
    dr = w.dr
    rl = w.rl
    r6 = w.r6
    dm = w.dm
    dm0 = w.dm0
    delta = w.delta
    fluxr = w.fluxr
    fluxu = w.fluxu
    fluxv = w.fluxv
    fluxw = w.fluxw
    fluxe = w.fluxe
    fluxq = w.fluxq

    paraset!(para)
    parabola!(nmin - 1, nmax + 1, para, sweeps.r, dr, r6, rl, sweeps.flat, w)
    parabola!(nmin - 1, nmax + 1, para, sweeps.u, du, u6, ul, sweeps.flat, w)
    parabola!(nmin - 1, nmax + 1, para, sweeps.v, dv, v6, vl, sweeps.flat, w)
    parabola!(nmin - 1, nmax + 1, para, sweeps.w, dw, w6, wl, sweeps.flat, w)
    parabola!(nmin - 1, nmax + 1, para, sweeps.q, dq, q6, ql, sweeps.flat, w)
    parabola!(nmin - 1, nmax + 1, para, sweeps.e, de, e6, el, sweeps.flat, w)

    for n in nmin:(nmax + 1)
        delta[n] = sweeps.xa[n] - sweeps.xa0[n]
    end

    for n in nmin:(nmax + 1)
        deltx = sweeps.xa[n] - sweeps.xa0[n]
        if deltx >= 0.0
            nn = n - 1
            fractn = 0.5 * deltx / sweeps.dx[nn]
            fractn2 = 1.0 - fourthd * fractn
            fluxr[n] = (rl[nn] + dr[nn] - fractn * (dr[nn] - fractn2 * r6[nn])) * delta[n]
            fluxu[n] = (ul[nn] + du[nn] - fractn * (du[nn] - fractn2 * u6[nn])) * fluxr[n]
            fluxv[n] = (vl[nn] + dv[nn] - fractn * (dv[nn] - fractn2 * v6[nn])) * fluxr[n]
            fluxw[n] = (wl[nn] + dw[nn] - fractn * (dw[nn] - fractn2 * w6[nn])) * fluxr[n]
            fluxe[n] = (el[nn] + de[nn] - fractn * (de[nn] - fractn2 * e6[nn])) * fluxr[n]
            fluxq[n] = (ql[nn] + dq[nn] - fractn * (dq[nn] - fractn2 * q6[nn])) * fluxr[n]
        else
            fractn = 0.5 * deltx / sweeps.dx[n]
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
        dm[n] = sweeps.r[n] * sweeps.dvol[n]
        dm0[n] = dm[n] + fluxr[n] - fluxr[n + 1]
        sweeps.r[n] = dm0[n] / sweeps.dvol0[n]
        sweeps.r[n] = max(smallr, sweeps.r[n])
        dm0[n] = 1.0 / (sweeps.r[n] * sweeps.dvol0[n])
        sweeps.u[n] = (sweeps.u[n] * dm[n] + fluxu[n] - fluxu[n + 1]) * dm0[n]
        sweeps.v[n] = (sweeps.v[n] * dm[n] + fluxv[n] - fluxv[n + 1]) * dm0[n]
        sweeps.w[n] = (sweeps.w[n] * dm[n] + fluxw[n] - fluxw[n + 1]) * dm0[n]
        sweeps.e[n] = (sweeps.e[n] * dm[n] + fluxe[n] - fluxe[n + 1]) * dm0[n]
        sweeps.q[n] = (sweeps.q[n] * dm[n] + fluxq[n] - fluxq[n + 1]) * dm0[n]

        ekin = 0.5 * (sweeps.u[n]^2 + sweeps.v[n]^2 + sweeps.w[n]^2)
        if ekin / sweeps.q[n] < 100.0
            sweeps.q[n] = sweeps.e[n] - ekin
        end

        sweeps.p[n] = gamm * sweeps.r[n] * sweeps.q[n]
        sweeps.p[n] = max(smallp, sweeps.p[n])
    end

    return nothing
end
