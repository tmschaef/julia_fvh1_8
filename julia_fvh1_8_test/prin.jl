# Julia port of prin.f90
# Writes 1D cuts of density, velocity, pressure, and entropy to output files.

include("vh1mods.jl")
include("zonemod.jl")
include("stress.jl")

using .Global
using .Zone

function myprin!(io_units::Dict{Int,<:IO})
    center_j = Zone.jmax ÷ 2 + 1
    center_k = Zone.kmax ÷ 2 + 1
    center_i = Zone.imax ÷ 2 + 1

    for i in 1:Zone.imax
        println(io_units[65], Global.time, " ", Zone.zxa[i], " ", Zone.zro[i, center_j, center_k])
    end

    for j in 1:Zone.jmax
        println(io_units[66], Global.time, " ", Zone.zya[j], " ", Zone.zro[center_i, j, center_k])
    end

    for k in 1:Zone.kmax
        println(io_units[67], Global.time, " ", Zone.zza[k], " ", Zone.zro[center_i, center_j, k])
    end

    for i in 1:Zone.imax
        println(io_units[68], Global.time, " ", Zone.zxa[i], " ", Zone.zux[i, center_j, center_k])
    end

    for j in 1:Zone.jmax
        println(io_units[88], Global.time, " ", Zone.zya[j], " ", Zone.zuy[center_i, j, center_k])
    end

    for k in 1:Zone.kmax
        println(io_units[98], Global.time, " ", Zone.zza[k], " ", Zone.zuz[center_i, center_j, k])
    end

    for i in 1:Zone.imax
        println(io_units[69], Global.time, " ", Zone.zxa[i], " ", Zone.zpr[i, center_j, center_k] / Zone.zro[i, center_j, center_k])
    end

    for i in 1:Zone.imax
        println(io_units[70], Global.time, " ", Zone.zxa[i], " ", Zone.zpr[i, center_j, center_k], " ", Zone.dpi[1, 1, i, center_j, center_k])
    end

    for i in 1:Zone.imax
        rloc = Zone.zro[i, center_j, center_k]
        ploc = Zone.zpr[i, center_j, center_k]
        sloc = slocal(rloc, ploc)
        println(io_units[71], Global.time, " ", Zone.zxa[i], " ", sloc, " ", sloc / rloc)
    end

    return nothing
end
