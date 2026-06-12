# Julia port of sweepx.f90
# Performs 1D X-direction hydro sweeps using the existing julia_port modules.

include("ppmlr.jl")
include("zonemod.jl")

using .Global
using .Zone
using .Sweeps

function sweepx!()
    Sweeps.nleft = Zone.nleftx
    Sweeps.nright = Zone.nrightx
    Sweeps.nmin = 7
    Sweeps.nmax = Zone.imax + 6
    
    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                n = i + 6

                Sweeps.r[n] = Zone.zro[i, j, k]
                Sweeps.p[n] = Zone.zpr[i, j, k]
                Sweeps.u[n] = Zone.zux[i, j, k]
                Sweeps.v[n] = Zone.zuy[i, j, k]
                Sweeps.w[n] = Zone.zuz[i, j, k]

                Sweeps.xa0[n] = Zone.zxa[i]
                Sweeps.dx0[n] = Zone.zdx[i]
                Sweeps.xa[n] = Zone.zxa[i]
                Sweeps.dx[n] = Zone.zdx[i]

                Sweeps.p[n] = max(Global.smallp, Sweeps.p[n])
                Sweeps.e[n] = Sweeps.p[n] / (Sweeps.r[n] * Global.gamm) + 0.5 * (Sweeps.u[n]^2 + Sweeps.v[n]^2 + Sweeps.w[n]^2)
            end
        
            ppmlr!()
        
            for i in 1:Zone.imax
                n = i + 6
                Zone.zro[i, j, k] = Sweeps.r[n]
                Zone.zpr[i, j, k] = Sweeps.p[n]
                Zone.zux[i, j, k] = Sweeps.u[n]
                Zone.zuy[i, j, k] = Sweeps.v[n]
                Zone.zuz[i, j, k] = Sweeps.w[n]
            end
        end
    end

    return nothing
end
