# Julia port of boundary.f90
# Applies the same ghost-cell boundary conditions as the Fortran routine.

include("vh1mods.jl")
include("zonemod.jl")

using .Global
using .Sweeps

"""
    boundary!()

Apply left and right boundary conditions to the 1D sweep arrays.
"""
function boundary!()
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    nleft = Sweeps.nleft
    nright = Sweeps.nright

    # Boundary-condition flags:
    #   nleft, nright = 0 : reflecting
    #                  = 1 : outflow (zero gradients)
    #                  = 2 : fixed inflow
    #                  = 3 : periodic

    for n in 1:6
        if nleft == 0
            Sweeps.dx[nmin - n] = Sweeps.dx[nmin + n - 1]
            Sweeps.xa[nmin - n] = Sweeps.xa[nmin - n + 1] - Sweeps.dx[nmin - n]
            Sweeps.dx0[nmin - n] = Sweeps.dx0[nmin + n - 1]
            Sweeps.xa0[nmin - n] = Sweeps.xa0[nmin - n + 1] - Sweeps.dx0[nmin - n]
            Sweeps.r[nmin - n] = Sweeps.r[nmin + n - 1]
            Sweeps.u[nmin - n] = -Sweeps.u[nmin + n - 1]
            Sweeps.v[nmin - n] = Sweeps.v[nmin + n - 1]
            Sweeps.w[nmin - n] = Sweeps.w[nmin + n - 1]
            Sweeps.p[nmin - n] = Sweeps.p[nmin + n - 1]
            Sweeps.e[nmin - n] = Sweeps.e[nmin + n - 1]
        elseif nleft == 1
            Sweeps.dx[nmin - n] = Sweeps.dx[nmin]
            Sweeps.xa[nmin - n] = Sweeps.xa[nmin - n + 1] - Sweeps.dx[nmin - n]
            Sweeps.dx0[nmin - n] = Sweeps.dx0[nmin]
            Sweeps.xa0[nmin - n] = Sweeps.xa0[nmin - n + 1] - Sweeps.dx0[nmin - n]
            Sweeps.r[nmin - n] = Sweeps.r[nmin]
            Sweeps.u[nmin - n] = Sweeps.u[nmin]
            Sweeps.v[nmin - n] = Sweeps.v[nmin]
            Sweeps.w[nmin - n] = Sweeps.w[nmin]
            Sweeps.p[nmin - n] = Sweeps.p[nmin]
            Sweeps.e[nmin - n] = Sweeps.e[nmin]
        elseif nleft == 2
            Sweeps.dx[nmin - n] = Sweeps.dx[nmin]
            Sweeps.xa[nmin - n] = Sweeps.xa[nmin - n + 1] - Sweeps.dx[nmin - n]
            Sweeps.dx0[nmin - n] = Sweeps.dx0[nmin]
            Sweeps.xa0[nmin - n] = Sweeps.xa0[nmin - n + 1] - Sweeps.dx0[nmin - n]
            Sweeps.r[nmin - n] = Global.dinflo
            Sweeps.u[nmin - n] = Global.uinflo
            Sweeps.v[nmin - n] = Global.vinflo
            Sweeps.w[nmin - n] = Global.winflo
            Sweeps.p[nmin - n] = Global.pinflo
            Sweeps.e[nmin - n] = Global.pinflo / (Global.dinflo * Global.gamm) + 0.5 * (Global.uinflo^2 + Global.vinflo^2 + Global.winflo^2)
        elseif nleft == 3
            Sweeps.dx[nmin - n] = Sweeps.dx[nmax + 1 - n]
            Sweeps.xa[nmin - n] = Sweeps.xa[nmin - n + 1] - Sweeps.dx[nmin - n]
            Sweeps.dx0[nmin - n] = Sweeps.dx0[nmax + 1 - n]
            Sweeps.xa0[nmin - n] = Sweeps.xa0[nmin - n + 1] - Sweeps.dx0[nmin - n]
            Sweeps.r[nmin - n] = Sweeps.r[nmax + 1 - n]
            Sweeps.u[nmin - n] = Sweeps.u[nmax + 1 - n]
            Sweeps.v[nmin - n] = Sweeps.v[nmax + 1 - n]
            Sweeps.w[nmin - n] = Sweeps.w[nmax + 1 - n]
            Sweeps.p[nmin - n] = Sweeps.p[nmax + 1 - n]
            Sweeps.e[nmin - n] = Sweeps.e[nmax + 1 - n]
        end
    end

    for n in 1:6
        if nright == 0
            Sweeps.dx[nmax + n] = Sweeps.dx[nmax + 1 - n]
            Sweeps.xa[nmax + n] = Sweeps.xa[nmax + n - 1] + Sweeps.dx[nmax + n - 1]
            Sweeps.dx0[nmax + n] = Sweeps.dx0[nmax + 1 - n]
            Sweeps.xa0[nmax + n] = Sweeps.xa0[nmax + n - 1] + Sweeps.dx0[nmax + n - 1]
            Sweeps.r[nmax + n] = Sweeps.r[nmax + 1 - n]
            Sweeps.u[nmax + n] = -Sweeps.u[nmax + 1 - n]
            Sweeps.v[nmax + n] = Sweeps.v[nmax + 1 - n]
            Sweeps.w[nmax + n] = Sweeps.w[nmax + 1 - n]
            Sweeps.p[nmax + n] = Sweeps.p[nmax + 1 - n]
            Sweeps.e[nmax + n] = Sweeps.e[nmax + 1 - n]
        elseif nright == 1
            Sweeps.dx[nmax + n] = Sweeps.dx[nmax]
            Sweeps.xa[nmax + n] = Sweeps.xa[nmax + n - 1] + Sweeps.dx[nmax + n - 1]
            Sweeps.dx0[nmax + n] = Sweeps.dx0[nmax]
            Sweeps.xa0[nmax + n] = Sweeps.xa0[nmax + n - 1] + Sweeps.dx0[nmax + n - 1]
            Sweeps.r[nmax + n] = Sweeps.r[nmax]
            Sweeps.u[nmax + n] = Sweeps.u[nmax]
            Sweeps.v[nmax + n] = Sweeps.v[nmax]
            Sweeps.w[nmax + n] = Sweeps.w[nmax]
            Sweeps.p[nmax + n] = Sweeps.p[nmax]
            Sweeps.e[nmax + n] = Sweeps.e[nmax]
        elseif nright == 2
            Sweeps.dx[nmax + n] = Sweeps.dx[nmax]
            Sweeps.xa[nmax + n] = Sweeps.xa[nmax + n - 1] + Sweeps.dx[nmax + n - 1]
            Sweeps.dx0[nmax + n] = Sweeps.dx0[nmax]
            Sweeps.xa0[nmax + n] = Sweeps.xa0[nmax + n - 1] + Sweeps.dx0[nmax + n - 1]
            Sweeps.r[nmax + n] = Global.dotflo
            Sweeps.u[nmax + n] = Global.uotflo
            Sweeps.v[nmax + n] = Global.votflo
            Sweeps.w[nmax + n] = Global.wotflo
            Sweeps.p[nmax + n] = Global.potflo
            Sweeps.e[nmax + n] = Global.potflo / (Global.dotflo * Global.gamm) + 0.5 * (Global.uotflo^2 + Global.votflo^2 + Global.wotflo^2)
        elseif nright == 3
            Sweeps.dx[nmax + n] = Sweeps.dx[nmin + n - 1]
            Sweeps.xa[nmax + n] = Sweeps.xa[nmax + n - 1] + Sweeps.dx[nmax + n - 1]
            Sweeps.dx0[nmax + n] = Sweeps.dx0[nmin + n - 1]
            Sweeps.xa0[nmax + n] = Sweeps.xa0[nmax + n - 1] + Sweeps.dx0[nmax + n - 1]
            Sweeps.r[nmax + n] = Sweeps.r[nmin + n - 1]
            Sweeps.u[nmax + n] = Sweeps.u[nmin + n - 1]
            Sweeps.v[nmax + n] = Sweeps.v[nmin + n - 1]
            Sweeps.w[nmax + n] = Sweeps.w[nmin + n - 1]
            Sweeps.p[nmax + n] = Sweeps.p[nmin + n - 1]
            Sweeps.e[nmax + n] = Sweeps.e[nmin + n - 1]
        end
    end

    return nothing
end
