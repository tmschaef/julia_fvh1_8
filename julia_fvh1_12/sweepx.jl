# Julia port of sweepx.f90 — parallel over (j, k) lines via Threads.@threads

function sweepx!(pool::SweepThreadPool, g::globalvars, z::zone, zg::zonegrid)
    nleft = zg.nleftx
    nright = zg.nrightx
    nmin = 7
    nmax = imax + 6

    Threads.@threads for jj in 1:(jmax * kmax)
        j = (jj - 1) ÷ kmax + 1
        k = (jj - 1) % kmax + 1
        tid = Threads.threadid()
        sweeps = pool.sweeps[tid]
        w = pool.work[tid]

        sweeps.nleft = nleft
        sweeps.nright = nright
        sweeps.nmin = nmin
        sweeps.nmax = nmax

        @inbounds for i in 1:imax
            n = i + 6

            sweeps.r[n] = z.zro[i, j, k]
            sweeps.p[n] = z.zpr[i, j, k]
            sweeps.u[n] = z.zux[i, j, k]
            sweeps.v[n] = z.zuy[i, j, k]
            sweeps.w[n] = z.zuz[i, j, k]

            sweeps.xa0[n] = zg.zxa[i]
            sweeps.dx0[n] = zg.zdx[i]
            sweeps.xa[n] = zg.zxa[i]
            sweeps.dx[n] = zg.zdx[i]

            sweeps.p[n] = max(smallp, sweeps.p[n])
            sweeps.e[n] = sweeps.p[n] / (sweeps.r[n] * gamm) +
                0.5 * (sweeps.u[n]^2 + sweeps.v[n]^2 + sweeps.w[n]^2)
        end

        ppmlr!(sweeps, g, w)

        @inbounds for i in 1:imax
            n = i + 6
            z.zro[i, j, k] = sweeps.r[n]
            z.zpr[i, j, k] = sweeps.p[n]
            z.zux[i, j, k] = sweeps.u[n]
            z.zuy[i, j, k] = sweeps.v[n]
            z.zuz[i, j, k] = sweeps.w[n]
        end
    end

    return nothing
end
