# --------------------------------------------------------------------------
#
#    11111111      11      11    11      11                1111  
#  1111    1111    11      11    11      11              111111  
#  1111            11      11    11      11                1111  
#    11111111      11      11    1111111111    ==          1111  
#          1111     11    11     11      11                1111  
#          1111    11    11      11      11                1111  
#  1111    1111      11  11      11      11                1111  
#    11111111         1111       11      11            1111111111
#
#
#                        VIRGINIA HYDRODYNAMICS #1
#
# --------------------------------------------------------------------------
# Julia port of vhone.f90
# Main VH-1 driver loop for the julia_fvh1_12 code (parallel sweeps + mc3d).
cd(@__DIR__)

include("vh1state.jl")
include("parabola.jl")
include("boundary.jl")
include("states.jl")
include("evolve.jl")
include("remap.jl")
include("ppmlr.jl")
include("init.jl")
include("stress.jl")
include("prin.jl")
include("dtcon.jl")
include("mc3d.jl")
include("sweepx.jl")
include("sweepy.jl")
include("sweepz.jl")

using Dates
using Random

const g = globalvars()
const z = zone()
const zg = zonegrid()
const zo = zoneobs()
const sweep_pool = SweepThreadPool()
const mc_pool = McThreadPool()

function parse_indat()
    prefix_val = ""
    endtime_val = 0.0
    tprin_val = 0.0

    if isfile("indat")
        for line in eachline("indat")
            line = strip(line)
            if isempty(line) || startswith(line, "!") || startswith(line, "&") || startswith(line, "/")
                continue
            end
            for part in split(line, ',')
                kv = split(part, '=')
                if length(kv) == 2
                    key = strip(lowercase(kv[1]))
                    value = strip(kv[2])
                    value = strip(replace(value, r"/\s*$" => ""))
                    if key == "prefix"
                        prefix_val = strip(value, ['\'', '"'])
                    elseif key == "endtime"
                        endtime_val = parse(Float64, value)
                    elseif key == "tprin"
                        tprin_val = parse(Float64, value)
                    end
                end
            end
        end
    end

    return prefix_val, endtime_val, tprin_val
end

function open_output_files()
    return Dict{Int,IO}(
        65 => open("rho_x.dat", "w"),
        66 => open("rho_y.dat", "w"),
        67 => open("rho_z.dat", "w"),
        68 => open("vel_x.dat", "w"),
        88 => open("vel_y.dat", "w"),
        98 => open("vel_z.dat", "w"),
        69 => open("temp_tran.dat", "w"),
        70 => open("pixx_tran.dat", "w"),
        71 => open("s_x.dat", "w"),
        75 => open("energy.dat", "w"),
        76 => open("radii.dat", "w"),
        78 => open("entropy.dat", "w"),
        79 => open("momentum.dat", "w"),
        80 => open("pipi.dat", "w"),
        81 => open("pipiim.dat", "w"),
        82 => open("pi-spec.dat", "w"),
        83 => open("rhokx.dat", "w")
    )
end

function close_output_files(io_units::Dict{Int,<:IO})
    for io in values(io_units)
        close(io)
    end
    return nothing
end

function format_write(io::IO, args...)
    println(io, join(args, " "))
    return nothing
end

function vhone!()
    Random.seed!(123456)
    g.ntrial = 0
    g.nacc = 0

    prefix, g.endtime, g.tprin = parse_indat()

    if max(imax, jmax, kmax) + 12 > maxsweep
        error("maxsweep too small")
    end
    if ndim != 3
        error("This is a 3d code")
    end

    io_units = open_output_files()
    history_name = string(prefix, "hst")
    history_io = open(history_name, "w")
    println(history_io, "History File for VH-1 simulation run on ", Dates.format(now(), "mm/dd/yyyy"))
    println(history_io)

    init!(g, z, zg)
    myprin!(io_units, g, z, zg)

    nprin = floor(Int, (g.tprin + 0.001) / g.dt)
    icor = 40
    iequ = 1
    nk = 4
    iconfmax = 5000

    if icor > Int(g.endtime / 2.0 / g.dt) + iequ
        error("icor too large")
    end

    zo.spec .= 0.0
    zo.spec2 .= 0.0

    iconf = 0
    ispec = 0
    ncycend = 100000

    Arepixky = zeros(Float64, iconfmax, nk)
    Aimpixky = zeros(Float64, iconfmax, nk)
    cor = zeros(Float64, icor + 1, nk)
    cori = zeros(Float64, icor + 1, nk)
    cor2 = zeros(Float64, icor + 1, nk)
    cori2 = zeros(Float64, icor + 1, nk)

    while g.ncycle < ncycend
        iconf += 1
        g.ncycle += 2
        g.ncycp += 2
        println("conf = ", iconf, " t = ", g.time, " dt = ", g.dt)

        etot, ekin, eint = energy!(z, zg, 0.0, 0.0, 0.0)
        stot, rtot = entropy!(z, zg, 0.0, 0.0)
        rx, ry, rz = scalesize!(z, zg, 0.0, 0.0, 0.0)

        ptot = zeros(Float64, 3)
        pdipole = zeros(Float64, 3, 3)
        mass, ptot, pdipole = momentum!(z, zg, 0.0, ptot, pdipole)
        pixky!(z, zo)

        ispec += 1
        for j in 1:jmax
            zo.spec[j] += zo.repixky[j]^2 + zo.impixky[j]^2
            zo.spec2[j] += (zo.repixky[j]^2 + zo.impixky[j]^2)^2
        end
        for k in 1:nk
            Arepixky[iconf, k] = zo.repixky[k]
            Aimpixky[iconf, k] = zo.impixky[k]
        end

        format_write(io_units[75], g.time, etot, ekin, eint)
        format_write(io_units[76], g.time, rx, ry, rz)
        format_write(io_units[78], g.time, stot, rtot, stot / rtot)
        format_write(io_units[79], g.time, mass, ptot[1], pdipole[1, 2] / mass, zo.impixky[1])
        format_write(io_units[83], g.time, zo.rhokx[1], zo.rhokx[2], zo.rhokx[3])

        sweepx!(sweep_pool, g, z, zg)
        sweepy!(sweep_pool, g, z, zg)
        sweepz!(sweep_pool, g, z, zg)

        g.time += g.dt
        g.timep += g.dt

        sweepz!(sweep_pool, g, z, zg)
        sweepy!(sweep_pool, g, z, zg)
        sweepx!(sweep_pool, g, z, zg)

        g.time += g.dt
        g.timep += g.dt

        vcon!(g, z, zg)
        g.dt = g.dtfix

        mc3d!(g, z, mc_pool)
        mc3d!(g, z, mc_pool)

        vcon!(g, z, zg)

        if g.ncycp >= nprin
            myprin!(io_units, g, z, zg)
            g.ncycp = 0
        end

        if g.time > g.endtime
            g.ncycle = ncycend
        end
    end

    isample = 0
    for i in iequ:(iconf - icor)
        isample += 1
        for j in 0:icor
            for k in 1:nk
                recor = Arepixky[i, k] * Arepixky[i + j, k] + Aimpixky[i, k] * Aimpixky[i + j, k]
                imcor = Arepixky[i, k] * Aimpixky[i + j, k] - Aimpixky[i, k] * Arepixky[i + j, k]
                cor[j + 1, k] += recor
                cori[j + 1, k] += imcor
                cor2[j + 1, k] += recor^2
                cori2[j + 1, k] += imcor^2
            end
        end
    end

    for j in 0:icor
        for k in 1:nk
            cor[j+1,k]  = cor[j+1,k]/float(isample)
            cori[j+1,k] = cori[j+1,k]/float(isample)
            cor2[j+1,k] = cor2[j+1,k]/float(isample)
            cori2[j+1,k] = cori2[j+1,k]/float(isample)
        end
        cork1 = cor[j + 1, 1] / cor[1, 1]
        cork2 = cor[j + 1, 2] / cor[1, 2]
        cork3 = cor[j + 1, 3] / cor[1, 3]
        dc1 = sqrt((cor2[j + 1, 1] - cor[j + 1, 1]^2) / isample) / cor[1, 1]
        dc2 = sqrt((cor2[j + 1, 2] - cor[j + 1, 2]^2) / isample) / cor[1, 2]
        dc3 = sqrt((cor2[j + 1, 3] - cor[j + 1, 3]^2) / isample) / cor[1, 3]
        format_write(io_units[80], 2.0 * g.dt * j, cork1, dc1, cork2, dc2, cork3, dc3)

        cork1 = cori[j + 1, 1] / cor[1, 1]
        cork2 = cori[j + 1, 2] / cor[1, 2]
        cork3 = cori[j + 1, 3] / cor[1, 3]
        dc1 = sqrt((cori2[j + 1, 1] - cori[j + 1, 1]^2) / isample) / cor[1, 1]
        dc2 = sqrt((cori2[j + 1, 2] - cori[j + 1, 2]^2) / isample) / cor[1, 2]
        dc3 = sqrt((cori2[j + 1, 3] - cori[j + 1, 3]^2) / isample) / cor[1, 3]
        format_write(io_units[81], 2.0 * g.dt * j, cork1, dc1, cork2, dc2, cork3, dc3)
    end

    for j in 1:jmax
        zo.spec[j] /= ispec
        zo.spec2[j] /= ispec
        format_write(io_units[82], j, zo.spec[j], sqrt((zo.spec2[j] - zo.spec[j]^2) / ispec))
    end

    println("acceptance rate ", Float64(g.nacc) / Float64(g.ntrial))

    close(history_io)
    close_output_files(io_units)

    return nothing
end
