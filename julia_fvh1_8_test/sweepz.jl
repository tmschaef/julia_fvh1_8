# Julia port of sweepz.f90
# Performs 1D Z-direction hydro sweeps using the existing julia_port modules.

include("ppmlr.jl")
include("zonemod.jl")

using .Global
using .Zone
using .Sweeps

function sweepz!()
    Sweeps.nleft = Zone.nleftz
    Sweeps.nright = Zone.nrightz
    Sweeps.nmin = 7
    Sweeps.nmax = Zone.kmax + 6

    for j in 1:Zone.jmax
        for i in 1:Zone.imax
            for k in 1:Zone.kmax
                n = k + 6

                Sweeps.r[n] = Zone.zro[i, j, k]
                Sweeps.p[n] = Zone.zpr[i, j, k]
                Sweeps.u[n] = Zone.zuz[i, j, k]
                Sweeps.v[n] = Zone.zux[i, j, k]
                Sweeps.w[n] = Zone.zuy[i, j, k]

                Sweeps.xa0[n] = Zone.zza[k]
                Sweeps.dx0[n] = Zone.zdz[k]
                Sweeps.xa[n] = Zone.zza[k]
                Sweeps.dx[n] = Zone.zdz[k]

                Sweeps.p[n] = max(Global.smallp, Sweeps.p[n])
                Sweeps.e[n] = Sweeps.p[n] / (Sweeps.r[n] * Global.gamm) + 0.5 * (Sweeps.u[n]^2 + Sweeps.v[n]^2 + Sweeps.w[n]^2)
            end

            ppmlr!()

            for k in 1:Zone.kmax
                n = k + 6
                Zone.zro[i, j, k] = Sweeps.r[n]
                Zone.zpr[i, j, k] = Sweeps.p[n]
                Zone.zuz[i, j, k] = Sweeps.u[n]
                Zone.zux[i, j, k] = Sweeps.v[n]
                Zone.zuy[i, j, k] = Sweeps.w[n]
            end
        end
    end

    return nothing
end
