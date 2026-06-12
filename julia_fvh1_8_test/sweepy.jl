# Julia port of sweepy.f90
# Performs 1D Y-direction hydro sweeps using the existing julia_port modules.

include("ppmlr.jl")
include("zonemod.jl")

using .Global
using .Zone
using .Sweeps

function sweepy!()
    Sweeps.nleft = Zone.nlefty
    Sweeps.nright = Zone.nrighty
    Sweeps.nmin = 7
    Sweeps.nmax = Zone.jmax + 6

    for k in 1:Zone.kmax
        for i in 1:Zone.imax
            for j in 1:Zone.jmax
                n = j + 6

                Sweeps.r[n] = Zone.zro[i, j, k]
                Sweeps.p[n] = Zone.zpr[i, j, k]
                Sweeps.u[n] = Zone.zuy[i, j, k]
                Sweeps.v[n] = Zone.zuz[i, j, k]
                Sweeps.w[n] = Zone.zux[i, j, k]

                Sweeps.xa0[n] = Zone.zya[j]
                Sweeps.dx0[n] = Zone.zdy[j]
                Sweeps.xa[n] = Zone.zya[j]
                Sweeps.dx[n] = Zone.zdy[j]

                Sweeps.p[n] = max(Global.smallp, Sweeps.p[n])
                Sweeps.e[n] = Sweeps.p[n] / (Sweeps.r[n] * Global.gamm) + 0.5 * (Sweeps.u[n]^2 + Sweeps.v[n]^2 + Sweeps.w[n]^2)
            end

            ppmlr!()

            for j in 1:Zone.jmax
                n = j + 6
                Zone.zro[i, j, k] = Sweeps.r[n]
                Zone.zpr[i, j, k] = Sweeps.p[n]
                Zone.zuy[i, j, k] = Sweeps.u[n]
                Zone.zuz[i, j, k] = Sweeps.v[n]
                Zone.zux[i, j, k] = Sweeps.w[n]
            end
        end
    end

    return nothing
end
