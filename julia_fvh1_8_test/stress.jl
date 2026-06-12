# Julia port of stress.f90
# Uses the already-ported modules in vh1mods.jl and zonemod.jl.

include("vh1mods.jl")
include("zonemod.jl")

using .Global
using .Zone

function energy!(etot::Real, ekin::Real, eint::Real)
    ekin = 0.0
    eint = 0.0

    for i in 1:Zone.imax
        for j in 1:Zone.jmax
            for k in 1:Zone.kmax
                d3x = Zone.zdx[i] * Zone.zdy[j] * Zone.zdz[k]
                v2 = Zone.zux[i, j, k]^2 + Zone.zuy[i, j, k]^2 + Zone.zuz[i, j, k]^2
                ekin += 0.5 * Zone.zro[i, j, k] * v2 * d3x
                eint += Zone.zpr[i, j, k] / Global.gamm * d3x
            end
        end
    end

    etot = ekin + eint
    return etot, ekin, eint
end

function scalesize!(rx::Real, ry::Real, rz::Real)
    xmtot = 0.0
    rx2 = 0.0
    ry2 = 0.0
    rz2 = 0.0

    for i in 1:Zone.imax
        for j in 1:Zone.jmax
            for k in 1:Zone.kmax
                xmtot += Zone.zro[i, j, k]
                rx2 += Zone.zro[i, j, k] * Zone.zxa[i]^2
                ry2 += Zone.zro[i, j, k] * Zone.zya[j]^2
                rz2 += Zone.zro[i, j, k] * Zone.zza[k]^2
            end
        end
    end

    rx = sqrt(rx2 / xmtot) * sqrt(2.0)
    ry = sqrt(ry2 / xmtot) * sqrt(2.0)
    rz = sqrt(rz2 / xmtot) * sqrt(2.0)
    return rx, ry, rz
end

function entropy!(stot::Real, rtot::Real)
    const_ = 400.0
    stot_sum = 0.0
    rtot_sum = 0.0

    for i in 1:Zone.imax
        for j in 1:Zone.jmax
            for k in 1:Zone.kmax
                d3x = Zone.zdx[i] * Zone.zdy[j] * Zone.zdz[k]
                rloc = Zone.zro[i, j, k]
                ploc = Zone.zpr[i, j, k]
                sloc = rloc * (2.5 + log(const_ * ploc^1.5 / rloc^2.5))
                rtot_sum += rloc * d3x
                stot_sum += sloc * d3x
            end
        end
    end

    stot = stot_sum
    rtot = rtot_sum
    return stot, rtot
end

function slocal(rloc::Real, ploc::Real)
    const_ = 400.0
    return rloc * (2.5 + log(const_ * ploc^1.5 / rloc^2.5))
end

function momentum!(mass::Real, ptot::Vector{Float64}, pdipole::Matrix{Float64})
    mass = 0.0
    fill!(ptot, 0.0)
    fill!(pdipole, 0.0)

    for i in 1:Zone.imax
        for j in 1:Zone.jmax
            for k in 1:Zone.kmax
                d3x = Zone.zdx[i] * Zone.zdy[j] * Zone.zdz[k]
                mass += d3x * Zone.zro[i, j, k]
                ptot[1] += d3x * Zone.zro[i, j, k] * Zone.zux[i, j, k]
                ptot[2] += d3x * Zone.zro[i, j, k] * Zone.zuy[i, j, k]
                ptot[3] += d3x * Zone.zro[i, j, k] * Zone.zuz[i, j, k]
                pdipole[1, 1] += d3x * Zone.zro[i, j, k] * Zone.zux[i, j, k] * Zone.zxa[i]
                pdipole[1, 2] += d3x * Zone.zro[i, j, k] * Zone.zux[i, j, k] * Zone.zya[j]
                pdipole[1, 3] += d3x * Zone.zro[i, j, k] * Zone.zux[i, j, k] * Zone.zza[k]
                pdipole[2, 1] += d3x * Zone.zro[i, j, k] * Zone.zuy[i, j, k] * Zone.zxa[i]
                pdipole[2, 2] += d3x * Zone.zro[i, j, k] * Zone.zuy[i, j, k] * Zone.zya[j]
                pdipole[2, 3] += d3x * Zone.zro[i, j, k] * Zone.zuy[i, j, k] * Zone.zza[k]
                pdipole[3, 1] += d3x * Zone.zro[i, j, k] * Zone.zuz[i, j, k] * Zone.zxa[i]
                pdipole[3, 2] += d3x * Zone.zro[i, j, k] * Zone.zuz[i, j, k] * Zone.zya[j]
                pdipole[3, 3] += d3x * Zone.zro[i, j, k] * Zone.zuz[i, j, k] * Zone.zza[k]
            end
        end
    end

    return mass, ptot, pdipole
end

function pixky!()
    nk = Zone.jmax
    fill!(Zone.repixky, 0.0)
    fill!(Zone.impixky, 0.0)
    fill!(Zone.rhokx, 0.0)

    for i in 1:Zone.imax
        for j in 1:Zone.jmax
            for k in 1:Zone.kmax
                pix = Zone.zro[i, j, k] * Zone.zux[i, j, k]
                for in_ in 1:nk
                    Zone.repixky[in_] += pix * cos(2.0 * π * j * in_ / Zone.jmax)
                    Zone.impixky[in_] += pix * sin(2.0 * π * j * in_ / Zone.jmax)
                    Zone.rhokx[in_] += Zone.zro[i, j, k] * sin(2.0 * π * (i - 1) * in_ / Zone.imax)
                end
            end
        end
    end

    xnorm = sqrt(Float64(Zone.imax * Zone.jmax * Zone.kmax))
    for in_ in 1:nk
        Zone.repixky[in_] /= xnorm
        Zone.impixky[in_] /= xnorm
        Zone.rhokx[in_] /= xnorm
    end

    return nothing
end
