# Julia port of states.f90
# Computes left/right characteristic states for the Riemann solver.

include("vh1mods.jl")
include("zonemod.jl")

using .Global
using .Sweeps
using .Sweepsize

"""
    states!(pl, ul, rl, p6, u6, r6, dp, du, dr,
            plft, ulft, rlft, prgh, urgh, rrgh)

Compute characteristic left/right states from zone-edge variables and slopes.
"""
function states!(pl::AbstractVector{<:Real}, ul::AbstractVector{<:Real}, rl::AbstractVector{<:Real},
                 p6::AbstractVector{<:Real}, u6::AbstractVector{<:Real}, r6::AbstractVector{<:Real},
                 dp::AbstractVector{<:Real}, du::AbstractVector{<:Real}, dr::AbstractVector{<:Real},
                 plft::AbstractVector{<:Real}, ulft::AbstractVector{<:Real}, rlft::AbstractVector{<:Real},
                 prgh::AbstractVector{<:Real}, urgh::AbstractVector{<:Real}, rrgh::AbstractVector{<:Real})
    nmin = Sweeps.nmin
    nmax = Sweeps.nmax
    fourthd = 4.0 / 3.0
    hdt = 0.5 * Global.dt

    Cdtdx = zeros(Float64, Sweepsize.maxsweep)
    fCdtdx = zeros(Float64, Sweepsize.maxsweep)

    for n in (nmin - 4):(nmax + 4)
        Cdtdx[n] = sqrt(Global.gam * Sweeps.p[n] / Sweeps.r[n]) / (Sweeps.dx[n] * Sweeps.radius)
        Cdtdx[n] = Cdtdx[n] * hdt
        fCdtdx[n] = 1.0 - fourthd * Cdtdx[n]
    end

    for n in (nmin - 4):(nmax + 4)
        np = n + 1
        plft[np] = pl[n] + dp[n] - Cdtdx[n] * (dp[n] - fCdtdx[n] * p6[n])
        ulft[np] = ul[n] + du[n] - Cdtdx[n] * (du[n] - fCdtdx[n] * u6[n])
        rlft[np] = rl[n] + dr[n] - Cdtdx[n] * (dr[n] - fCdtdx[n] * r6[n])
        plft[np] = max(Global.smallp, plft[np])
        rlft[np] = max(Global.smallr, rlft[np])

        prgh[n] = pl[n] + Cdtdx[n] * (dp[n] + fCdtdx[n] * p6[n])
        urgh[n] = ul[n] + Cdtdx[n] * (du[n] + fCdtdx[n] * u6[n])
        rrgh[n] = rl[n] + Cdtdx[n] * (dr[n] + fCdtdx[n] * r6[n])
        prgh[n] = max(Global.smallp, prgh[n])
        rrgh[n] = max(Global.smallr, rrgh[n])
    end
    
    return nothing
end
