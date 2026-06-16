# Julia port of stress.f90

function energy!(z::zone, zg::zonegrid, etot::Real, ekin::Real, eint::Real)
    ekin = 0.0
    eint = 0.0

    for i in 1:imax
        for j in 1:jmax
            for k in 1:kmax
                d3x = zg.zdx[i] * zg.zdy[j] * zg.zdz[k]
                v2 = z.zux[i, j, k]^2 + z.zuy[i, j, k]^2 + z.zuz[i, j, k]^2
                ekin += 0.5 * z.zro[i, j, k] * v2 * d3x
                eint += z.zpr[i, j, k] / gamm * d3x
            end
        end
    end

    etot = ekin + eint
    return etot, ekin, eint
end

function scalesize!(z::zone, zg::zonegrid, rx::Real, ry::Real, rz::Real)
    xmtot = 0.0
    rx2 = 0.0
    ry2 = 0.0
    rz2 = 0.0

    for i in 1:imax
        for j in 1:jmax
            for k in 1:kmax
                xmtot += z.zro[i, j, k]
                rx2 += z.zro[i, j, k] * zg.zxa[i]^2
                ry2 += z.zro[i, j, k] * zg.zya[j]^2
                rz2 += z.zro[i, j, k] * zg.zza[k]^2
            end
        end
    end

    rx = sqrt(rx2 / xmtot) * sqrt(2.0)
    ry = sqrt(ry2 / xmtot) * sqrt(2.0)
    rz = sqrt(rz2 / xmtot) * sqrt(2.0)
    return rx, ry, rz
end

function entropy!(z::zone, zg::zonegrid, stot::Real, rtot::Real)
    const_ = 400.0
    stot_sum = 0.0
    rtot_sum = 0.0

    for i in 1:imax
        for j in 1:jmax
            for k in 1:kmax
                d3x = zg.zdx[i] * zg.zdy[j] * zg.zdz[k]
                rloc = z.zro[i, j, k]
                ploc = z.zpr[i, j, k]
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

function momentum!(z::zone, zg::zonegrid, mass::Real, ptot::Vector{Float64}, pdipole::Matrix{Float64})
    mass = 0.0
    fill!(ptot, 0.0)
    fill!(pdipole, 0.0)

    for i in 1:imax
        for j in 1:jmax
            for k in 1:kmax
                d3x = zg.zdx[i] * zg.zdy[j] * zg.zdz[k]
                mass += d3x * z.zro[i, j, k]
                ptot[1] += d3x * z.zro[i, j, k] * z.zux[i, j, k]
                ptot[2] += d3x * z.zro[i, j, k] * z.zuy[i, j, k]
                ptot[3] += d3x * z.zro[i, j, k] * z.zuz[i, j, k]
                pdipole[1, 1] += d3x * z.zro[i, j, k] * z.zux[i, j, k] * zg.zxa[i]
                pdipole[1, 2] += d3x * z.zro[i, j, k] * z.zux[i, j, k] * zg.zya[j]
                pdipole[1, 3] += d3x * z.zro[i, j, k] * z.zux[i, j, k] * zg.zza[k]
                pdipole[2, 1] += d3x * z.zro[i, j, k] * z.zuy[i, j, k] * zg.zxa[i]
                pdipole[2, 2] += d3x * z.zro[i, j, k] * z.zuy[i, j, k] * zg.zya[j]
                pdipole[2, 3] += d3x * z.zro[i, j, k] * z.zuy[i, j, k] * zg.zza[k]
                pdipole[3, 1] += d3x * z.zro[i, j, k] * z.zuz[i, j, k] * zg.zxa[i]
                pdipole[3, 2] += d3x * z.zro[i, j, k] * z.zuz[i, j, k] * zg.zya[j]
                pdipole[3, 3] += d3x * z.zro[i, j, k] * z.zuz[i, j, k] * zg.zza[k]
            end
        end
    end

    return mass, ptot, pdipole
end

function pixky!(z::zone, zo::zoneobs)
    nk = jmax

    Threads.@threads for in_ in 1:nk
        rep = 0.0
        imp = 0.0
        rh = 0.0
        @inbounds for k in 1:kmax, j in 1:jmax, i in 1:imax
            rho = z.zro[i, j, k]
            pix = rho * z.zux[i, j, k]
            rep += pix * PIXKY_COS_J[j, in_]
            imp += pix * PIXKY_SIN_J[j, in_]
            rh += rho * PIXKY_SIN_I[i, in_]
        end
        zo.repixky[in_] = rep * PIXKY_XNORM
        zo.impixky[in_] = imp * PIXKY_XNORM
        zo.rhokx[in_] = rh * PIXKY_XNORM
    end

    return nothing
end
