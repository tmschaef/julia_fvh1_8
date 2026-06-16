# Julia port of states.f90

"""
    states!(sweeps, g, pl, ul, rl, p6, u6, r6, dp, du, dr,
            plft, ulft, rlft, prgh, urgh, rrgh, w)

Compute characteristic left/right states from zone-edge variables and slopes.
"""
function states!(sweeps::Sweeps, g::globalvars,
                 pl::AbstractVector{<:Real}, ul::AbstractVector{<:Real}, rl::AbstractVector{<:Real},
                 p6::AbstractVector{<:Real}, u6::AbstractVector{<:Real}, r6::AbstractVector{<:Real},
                 dp::AbstractVector{<:Real}, du::AbstractVector{<:Real}, dr::AbstractVector{<:Real},
                 plft::AbstractVector{<:Real}, ulft::AbstractVector{<:Real}, rlft::AbstractVector{<:Real},
                 prgh::AbstractVector{<:Real}, urgh::AbstractVector{<:Real}, rrgh::AbstractVector{<:Real},
                 w::PpmWork)
    nmin = sweeps.nmin
    nmax = sweeps.nmax
    fourthd = 4.0 / 3.0
    hdt = 0.5 * g.dt

    Cdtdx = w.Cdtdx
    fCdtdx = w.fCdtdx

    for n in (nmin - 4):(nmax + 4)
        Cdtdx[n] = sqrt(gam * sweeps.p[n] / sweeps.r[n]) / (sweeps.dx[n] * sweeps.radius)
        Cdtdx[n] = Cdtdx[n] * hdt
        fCdtdx[n] = 1.0 - fourthd * Cdtdx[n]
    end

    for n in (nmin - 4):(nmax + 4)
        np = n + 1
        plft[np] = pl[n] + dp[n] - Cdtdx[n] * (dp[n] - fCdtdx[n] * p6[n])
        ulft[np] = ul[n] + du[n] - Cdtdx[n] * (du[n] - fCdtdx[n] * u6[n])
        rlft[np] = rl[n] + dr[n] - Cdtdx[n] * (dr[n] - fCdtdx[n] * r6[n])
        plft[np] = max(smallp, plft[np])
        rlft[np] = max(smallr, rlft[np])

        prgh[n] = pl[n] + Cdtdx[n] * (dp[n] + fCdtdx[n] * p6[n])
        urgh[n] = ul[n] + Cdtdx[n] * (du[n] + fCdtdx[n] * u6[n])
        rrgh[n] = rl[n] + Cdtdx[n] * (dr[n] + fCdtdx[n] * r6[n])
        prgh[n] = max(smallp, prgh[n])
        rrgh[n] = max(smallr, rrgh[n])
    end

    return nothing
end
