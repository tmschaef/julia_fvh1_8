# Julia port of sweepz.f90 — parallel over (i, j) lines via Threads.@threads

function sweepz!(pool::SweepThreadPool, g::globalvars, z::zone, zg::zonegrid)
    nleft = zg.nleftz
    nright = zg.nrightz
    nmin = 7
    nmax = kmax + 6

    Threads.@threads for jj in 1:(imax * jmax)
        i = (jj - 1) ÷ jmax + 1
        j = (jj - 1) % jmax + 1
        tid = Threads.threadid()
        sweeps = pool.sweeps[tid]
        w = pool.work[tid]

        sweeps.nleft = nleft
        sweeps.nright = nright
        sweeps.nmin = nmin
        sweeps.nmax = nmax

        @inbounds for k in 1:kmax
            n = k + 6

            sweeps.r[n] = z.zro[i, j, k]
            sweeps.p[n] = z.zpr[i, j, k]
            sweeps.u[n] = z.zuz[i, j, k]
            sweeps.v[n] = z.zux[i, j, k]
            sweeps.w[n] = z.zuy[i, j, k]

            sweeps.xa0[n] = zg.zza[k]
            sweeps.dx0[n] = zg.zdz[k]
            sweeps.xa[n] = zg.zza[k]
            sweeps.dx[n] = zg.zdz[k]

            sweeps.p[n] = max(smallp, sweeps.p[n])
            sweeps.e[n] = sweeps.p[n] / (sweeps.r[n] * gamm) +
                0.5 * (sweeps.u[n]^2 + sweeps.v[n]^2 + sweeps.w[n]^2)
        end

        ppmlr!(sweeps, g, w)

        @inbounds for k in 1:kmax
            n = k + 6
            z.zro[i, j, k] = sweeps.r[n]
            z.zpr[i, j, k] = sweeps.p[n]
            z.zuz[i, j, k] = sweeps.u[n]
            z.zux[i, j, k] = sweeps.v[n]
            z.zuy[i, j, k] = sweeps.w[n]
        end
    end

    return nothing
end
