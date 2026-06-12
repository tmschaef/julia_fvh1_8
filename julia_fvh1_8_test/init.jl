# Julia port of init.f90
# Uses the already-ported modules from vh1mods.jl and zonemod.jl.

include("vh1mods.jl")
include("zonemod.jl")

using .Global
using .Zone

function grid!(nzones::Integer, xmin::Real, xmax::Real, xa::Vector{Float64}, xc::Vector{Float64}, dx::Vector{Float64})
    dxfac = (xmax - xmin) / Float64(nzones)
    for n in 1:nzones
        xa[n] = xmin + (n - 1) * dxfac
        dx[n] = dxfac
        xc[n] = xa[n] + 0.5 * dx[n]
    end
    return nothing
end

function init!()
    # Set up geometry and boundary conditions of grid
    Global.pi = 2.0 * asin(1.0)
    
    Zone.nleftx = 3
    Zone.nrightx = 3
    Zone.nlefty = 3
    Zone.nrighty = 3
    Zone.nleftz = 3
    Zone.nrightz = 3

    xmin = -5.0
    xmax = 5.0
    ymin = -5.0
    ymax = 5.0
    zmin = -5.0
    zmax = 5.0

    # Set time and cycle counters
    Global.time = 0.0
    Global.timep = 0.0
    Global.ncycle = 0
    Global.ncycp = 0

    # Set up grid coordinates
    grid!(Zone.imax, xmin, xmax, Zone.zxa, Zone.zxc, Zone.zdx)
    grid!(Zone.jmax, ymin, ymax, Zone.zya, Zone.zyc, Zone.zdy)
    grid!(Zone.kmax, zmin, zmax, Zone.zza, Zone.zzc, Zone.zdz)
    
    # Set up parameters for the problem
    Global.dtfix = 0.05
    Global.A0 = Zone.zdx[1] * Zone.zdy[1]
    Global.V0 = Global.A0 * Zone.zdz[1]

    # Box
    Global.tav = 1.0
    Global.nav = 400.0
    Global.etaav = 0.1 * Global.nav 
    
    # Trap
    lambda = 1.0
    e0ef = 1.0
    n00 = 1.0

    Global.gam = 5.0 / 3.0
    Global.gamm = Global.gam - 1.0

    # density, temperature dependence of eta
    Global.alpha = 0.0
    Global.alphat = 0.0

    # initialize profile, gaussian distribution
    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                rr2 = Zone.zxa[i]^2 + Zone.zya[j]^2 + lambda^2 * Zone.zza[k]^2
                vconf = rr2 / 1.0
                rhoi = Global.nav * exp(-vconf)
                Zone.zro[i, j, k] = max(rhoi, Global.smallr)
                Zone.zpr[i, j, k] = Zone.zro[i, j, k] * Global.tav
                Zone.zux[i, j, k] = 0.0
                Zone.zuy[i, j, k] = 0.0
                Zone.zuz[i, j, k] = 0.0
            end
        end
    end

    # initialize profile, homogenous density
    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                Zone.zro[i, j, k] = Global.nav
                rr2 = Zone.zxa[i]^2 + Zone.zya[j]^2 + Zone.zza[k]^2
                Zone.zro[i, j, k] = Global.nav * (1.0 + 0.0 * exp(-rr2))
                Zone.zpr[i, j, k] = Zone.zro[i, j, k] * Global.tav
                Zone.zux[i, j, k] = 0.0
                Zone.zuy[i, j, k] = 0.0
                Zone.zuz[i, j, k] = 0.0
            end
        end
    end

    # initialize profile, shear wave
    v0shear = 0.5
    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                Zone.zro[i, j, k] = Global.nav
                rr2 = Zone.zxa[i]^2 + Zone.zya[j]^2 + Zone.zza[k]^2
                Zone.zpr[i, j, k] = Zone.zro[i, j, k] * Global.tav
                Zone.zux[i, j, k] = v0shear * sin(Global.pi * Zone.zya[j] / ymax)
                Zone.zuy[i, j, k] = 0.0
                Zone.zuz[i, j, k] = 0.0
            end
        end
    end

    # initialize profile, density wave
    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                Zone.zro[i, j, k] = Global.nav
                Zone.zro[i, j, k] = Global.nav * (1.0 + 0.1 * sin(Global.pi * Zone.zxa[i] / xmax))
                Zone.zpr[i, j, k] = Global.nav * Global.tav * (Zone.zro[i, j, k] / Global.nav)^Global.gam
                Zone.zux[i, j, k] = 0.0
                Zone.zuy[i, j, k] = 0.0
                Zone.zuz[i, j, k] = 0.0
            end
        end
    end

    Global.dt = Global.dtfix
    return nothing
end
