# Julia port of boundary.f90

"""
    boundary!(sweeps::Sweeps)

Apply left and right boundary conditions to the 1D sweep arrays.
"""
function boundary!(sweeps::Sweeps)
    nmin = sweeps.nmin
    nmax = sweeps.nmax
    nleft = sweeps.nleft
    nright = sweeps.nright

    for n in 1:6
        if nleft == 0
            sweeps.dx[nmin - n] = sweeps.dx[nmin + n - 1]
            sweeps.xa[nmin - n] = sweeps.xa[nmin - n + 1] - sweeps.dx[nmin - n]
            sweeps.dx0[nmin - n] = sweeps.dx0[nmin + n - 1]
            sweeps.xa0[nmin - n] = sweeps.xa0[nmin - n + 1] - sweeps.dx0[nmin - n]
            sweeps.r[nmin - n] = sweeps.r[nmin + n - 1]
            sweeps.u[nmin - n] = -sweeps.u[nmin + n - 1]
            sweeps.v[nmin - n] = sweeps.v[nmin + n - 1]
            sweeps.w[nmin - n] = sweeps.w[nmin + n - 1]
            sweeps.p[nmin - n] = sweeps.p[nmin + n - 1]
            sweeps.e[nmin - n] = sweeps.e[nmin + n - 1]
        elseif nleft == 1
            sweeps.dx[nmin - n] = sweeps.dx[nmin]
            sweeps.xa[nmin - n] = sweeps.xa[nmin - n + 1] - sweeps.dx[nmin - n]
            sweeps.dx0[nmin - n] = sweeps.dx0[nmin]
            sweeps.xa0[nmin - n] = sweeps.xa0[nmin - n + 1] - sweeps.dx0[nmin - n]
            sweeps.r[nmin - n] = sweeps.r[nmin]
            sweeps.u[nmin - n] = sweeps.u[nmin]
            sweeps.v[nmin - n] = sweeps.v[nmin]
            sweeps.w[nmin - n] = sweeps.w[nmin]
            sweeps.p[nmin - n] = sweeps.p[nmin]
            sweeps.e[nmin - n] = sweeps.e[nmin]
        elseif nleft == 2
            sweeps.dx[nmin - n] = sweeps.dx[nmin]
            sweeps.xa[nmin - n] = sweeps.xa[nmin - n + 1] - sweeps.dx[nmin - n]
            sweeps.dx0[nmin - n] = sweeps.dx0[nmin]
            sweeps.xa0[nmin - n] = sweeps.xa0[nmin - n + 1] - sweeps.dx0[nmin - n]
            sweeps.r[nmin - n] = dinflo
            sweeps.u[nmin - n] = uinflo
            sweeps.v[nmin - n] = vinflo
            sweeps.w[nmin - n] = winflo
            sweeps.p[nmin - n] = pinflo
            sweeps.e[nmin - n] = pinflo / (dinflo * gamm) + 0.5 * (uinflo^2 + vinflo^2 + winflo^2)
        elseif nleft == 3
            sweeps.dx[nmin - n] = sweeps.dx[nmax + 1 - n]
            sweeps.xa[nmin - n] = sweeps.xa[nmin - n + 1] - sweeps.dx[nmin - n]
            sweeps.dx0[nmin - n] = sweeps.dx0[nmax + 1 - n]
            sweeps.xa0[nmin - n] = sweeps.xa0[nmin - n + 1] - sweeps.dx0[nmin - n]
            sweeps.r[nmin - n] = sweeps.r[nmax + 1 - n]
            sweeps.u[nmin - n] = sweeps.u[nmax + 1 - n]
            sweeps.v[nmin - n] = sweeps.v[nmax + 1 - n]
            sweeps.w[nmin - n] = sweeps.w[nmax + 1 - n]
            sweeps.p[nmin - n] = sweeps.p[nmax + 1 - n]
            sweeps.e[nmin - n] = sweeps.e[nmax + 1 - n]
        end
    end

    for n in 1:6
        if nright == 0
            sweeps.dx[nmax + n] = sweeps.dx[nmax + 1 - n]
            sweeps.xa[nmax + n] = sweeps.xa[nmax + n - 1] + sweeps.dx[nmax + n - 1]
            sweeps.dx0[nmax + n] = sweeps.dx0[nmax + 1 - n]
            sweeps.xa0[nmax + n] = sweeps.xa0[nmax + n - 1] + sweeps.dx0[nmax + n - 1]
            sweeps.r[nmax + n] = sweeps.r[nmax + 1 - n]
            sweeps.u[nmax + n] = -sweeps.u[nmax + 1 - n]
            sweeps.v[nmax + n] = sweeps.v[nmax + 1 - n]
            sweeps.w[nmax + n] = sweeps.w[nmax + 1 - n]
            sweeps.p[nmax + n] = sweeps.p[nmax + 1 - n]
            sweeps.e[nmax + n] = sweeps.e[nmax + 1 - n]
        elseif nright == 1
            sweeps.dx[nmax + n] = sweeps.dx[nmax]
            sweeps.xa[nmax + n] = sweeps.xa[nmax + n - 1] + sweeps.dx[nmax + n - 1]
            sweeps.dx0[nmax + n] = sweeps.dx0[nmax]
            sweeps.xa0[nmax + n] = sweeps.xa0[nmax + n - 1] + sweeps.dx0[nmax + n - 1]
            sweeps.r[nmax + n] = sweeps.r[nmax]
            sweeps.u[nmax + n] = sweeps.u[nmax]
            sweeps.v[nmax + n] = sweeps.v[nmax]
            sweeps.w[nmax + n] = sweeps.w[nmax]
            sweeps.p[nmax + n] = sweeps.p[nmax]
            sweeps.e[nmax + n] = sweeps.e[nmax]
        elseif nright == 2
            sweeps.dx[nmax + n] = sweeps.dx[nmax]
            sweeps.xa[nmax + n] = sweeps.xa[nmax + n - 1] + sweeps.dx[nmax + n - 1]
            sweeps.dx0[nmax + n] = sweeps.dx0[nmax]
            sweeps.xa0[nmax + n] = sweeps.xa0[nmax + n - 1] + sweeps.dx0[nmax + n - 1]
            sweeps.r[nmax + n] = dotflo
            sweeps.u[nmax + n] = uotflo
            sweeps.v[nmax + n] = votflo
            sweeps.w[nmax + n] = wotflo
            sweeps.p[nmax + n] = potflo
            sweeps.e[nmax + n] = potflo / (dotflo * gamm) + 0.5 * (uotflo^2 + votflo^2 + wotflo^2)
        elseif nright == 3
            sweeps.dx[nmax + n] = sweeps.dx[nmin + n - 1]
            sweeps.xa[nmax + n] = sweeps.xa[nmax + n - 1] + sweeps.dx[nmax + n - 1]
            sweeps.dx0[nmax + n] = sweeps.dx0[nmin + n - 1]
            sweeps.xa0[nmax + n] = sweeps.xa0[nmax + n - 1] + sweeps.dx0[nmax + n - 1]
            sweeps.r[nmax + n] = sweeps.r[nmin + n - 1]
            sweeps.u[nmax + n] = sweeps.u[nmin + n - 1]
            sweeps.v[nmax + n] = sweeps.v[nmin + n - 1]
            sweeps.w[nmax + n] = sweeps.w[nmin + n - 1]
            sweeps.p[nmax + n] = sweeps.p[nmin + n - 1]
            sweeps.e[nmax + n] = sweeps.e[nmin + n - 1]
        end
    end

    return nothing
end
