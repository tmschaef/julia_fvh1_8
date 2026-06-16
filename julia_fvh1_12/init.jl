# Julia port of init.f90

function grid!(nzones::Integer, xmin::Real, xmax::Real, xa::Vector{Float64}, xc::Vector{Float64}, dx::Vector{Float64})
    dxfac = (xmax - xmin) / Float64(nzones)
    for n in 1:nzones
        xa[n] = xmin + (n - 1) * dxfac
        dx[n] = dxfac
        xc[n] = xa[n] + 0.5 * dx[n]
    end
    return nothing
end

function init!(g::globalvars, z::zone, zg::zonegrid)
    zg.nleftx = 3
    zg.nrightx = 3
    zg.nlefty = 3
    zg.nrighty = 3
    zg.nleftz = 3
    zg.nrightz = 3

    xmin = -5.0
    xmax = 5.0
    ymin = -5.0
    ymax = 5.0
    zmin = -5.0
    zmax = 5.0

    g.time = 0.0
    g.timep = 0.0
    g.ncycle = 0
    g.ncycp = 0

    grid!(imax, xmin, xmax, zg.zxa, zg.zxc, zg.zdx)
    grid!(jmax, ymin, ymax, zg.zya, zg.zyc, zg.zdy)
    grid!(kmax, zmin, zmax, zg.zza, zg.zzc, zg.zdz)

    g.dtfix = 0.05

    lambda = 1.0
    e0ef = 1.0
    n00 = 1.0

    for k in 1:kmax
        for j in 1:jmax
            for i in 1:imax
                rr2 = zg.zxa[i]^2 + zg.zya[j]^2 + lambda^2 * zg.zza[k]^2
                vconf = rr2 / 1.0
                rhoi = nav * exp(-vconf)
                z.zro[i, j, k] = max(rhoi, smallr)
                z.zpr[i, j, k] = z.zro[i, j, k] * tav
                z.zux[i, j, k] = 0.0
                z.zuy[i, j, k] = 0.0
                z.zuz[i, j, k] = 0.0
            end
        end
    end

    for k in 1:kmax
        for j in 1:jmax
            for i in 1:imax
                z.zro[i, j, k] = nav
                rr2 = zg.zxa[i]^2 + zg.zya[j]^2 + zg.zza[k]^2
                z.zro[i, j, k] = nav * (1.0 + 0.0 * exp(-rr2))
                z.zpr[i, j, k] = z.zro[i, j, k] * tav
                z.zux[i, j, k] = 0.0
                z.zuy[i, j, k] = 0.0
                z.zuz[i, j, k] = 0.0
            end
        end
    end

    v0shear = 0.5
    for k in 1:kmax
        for j in 1:jmax
            for i in 1:imax
               z.zro[i, j, k] = nav
               rr2 = zg.zxa[i]^2 + zg.zya[j]^2 + zg.zza[k]^2
               z.zpr[i, j, k] = z.zro[i, j, k] * tav
               z.zux[i, j, k] = v0shear * sin(pi * zg.zya[j] / ymax)
               z.zuy[i, j, k] = 0.0
               z.zuz[i, j, k] = 0.0
            end
        end
    end

    for k in 1:kmax
        for j in 1:jmax
            for i in 1:imax
                z.zro[i, j, k] = nav
                z.zro[i, j, k] = nav * (1.0 + 0.1 * sin(pi * zg.zxa[i] / xmax))
                z.zpr[i, j, k] = nav * tav * (z.zro[i, j, k] / nav)^gam
                z.zux[i, j, k] = 0.0
                z.zuy[i, j, k] = 0.0
                z.zuz[i, j, k] = 0.0
            end
        end
    end

    g.dt = g.dtfix
    return nothing
end
