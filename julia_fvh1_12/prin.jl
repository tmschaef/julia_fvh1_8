# Julia port of prin.f90

function myprin!(io_units::Dict{Int,<:IO}, g::globalvars, z::zone, zg::zonegrid)
    center_j = jmax ÷ 2 + 1
    center_k = kmax ÷ 2 + 1
    center_i = imax ÷ 2 + 1

    for i in 1:imax
        println(io_units[65], g.time, " ", zg.zxa[i], " ", z.zro[i, center_j, center_k])
    end

    for j in 1:jmax
        println(io_units[66], g.time, " ", zg.zya[j], " ", z.zro[center_i, j, center_k])
    end

    for k in 1:kmax
        println(io_units[67], g.time, " ", zg.zza[k], " ", z.zro[center_i, center_j, k])
    end

    for i in 1:imax
        println(io_units[68], g.time, " ", zg.zxa[i], " ", z.zux[i, center_j, center_k])
    end

    for j in 1:jmax
        println(io_units[88], g.time, " ", zg.zya[j], " ", z.zuy[center_i, j, center_k])
    end

    for k in 1:kmax
        println(io_units[98], g.time, " ", zg.zza[k], " ", z.zuz[center_i, center_j, k])
    end

    for i in 1:imax
        println(io_units[69], g.time, " ", zg.zxa[i], " ", z.zpr[i, center_j, center_k] / z.zro[i, center_j, center_k])
    end

    for i in 1:imax
        println(io_units[70], g.time, " ", zg.zxa[i], " ", z.zpr[i, center_j, center_k])
    end

    for i in 1:imax
        rloc = z.zro[i, center_j, center_k]
        ploc = z.zpr[i, center_j, center_k]
        sloc = slocal(rloc, ploc)
        println(io_units[71], g.time, " ", zg.zxa[i], " ", sloc, " ", sloc / rloc)
    end

    return nothing
end
