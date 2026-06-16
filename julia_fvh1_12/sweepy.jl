# Julia port of sweepy.f90 — parallel over (i, k) lines via Threads.@threads

function sweepy!(pool::SweepThreadPool, g::globalvars, z::zone, zg::zonegrid)
    nleft = zg.nlefty
    nright = zg.nrighty
    nmin = 7
    nmax = jmax + 6

    Threads.@threads for jj in 1:(imax * kmax)
        i = (jj - 1) ÷ kmax + 1
        k = (jj - 1) % kmax + 1
        tid = Threads.threadid()
        sweeps = pool.sweeps[tid]
        w = pool.work[tid]

        sweeps.nleft = nleft
        sweeps.nright = nright
        sweeps.nmin = nmin
        sweeps.nmax = nmax

        @inbounds for j in 1:jmax
            n = j + 6

            sweeps.r[n] = z.zro[i, j, k]
            sweeps.p[n] = z.zpr[i, j, k]
            sweeps.u[n] = z.zuy[i, j, k]
            sweeps.v[n] = z.zuz[i, j, k]
            sweeps.w[n] = z.zux[i, j, k]

            sweeps.xa0[n] = zg.zya[j]
            sweeps.dx0[n] = zg.zdy[j]
            sweeps.xa[n] = zg.zya[j]
            sweeps.dx[n] = zg.zdy[j]

            sweeps.p[n] = max(smallp, sweeps.p[n])
            sweeps.e[n] = sweeps.p[n] / (sweeps.r[n] * gamm) +
                0.5 * (sweeps.u[n]^2 + sweeps.v[n]^2 + sweeps.w[n]^2)
        end

        ppmlr!(sweeps, g, w)

        @inbounds for j in 1:jmax
            n = j + 6
            z.zro[i, j, k] = sweeps.r[n]
            z.zpr[i, j, k] = sweeps.p[n]
            z.zuy[i, j, k] = sweeps.u[n]
            z.zuz[i, j, k] = sweeps.v[n]
            z.zux[i, j, k] = sweeps.w[n]
        end
    end

    return nothing
end
