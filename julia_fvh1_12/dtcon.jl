# Julia port of dtcon.f90

"""
    vcon!(g::globalvars, z::zone, zg::zonegrid)

Clamp the velocity components in the zone arrays to the CFL bounds implied by dt.
"""
function vcon!(g::globalvars, z::zone, zg::zonegrid)
    vxbound = zg.zdx[1] / g.dt * 0.45
    vybound = zg.zdy[1] / g.dt * 0.45
    vzbound = zg.zdz[1] / g.dt * 0.45

    for k in 1:kmax
        for j in 1:jmax
            for i in 1:imax
                xvel = abs(z.zux[i, j, k])
                yvel = abs(z.zuy[i, j, k])
                zvel = abs(z.zuz[i, j, k])

                if xvel > vxbound
                    z.zux[i, j, k] *= vxbound / xvel
                end
                if yvel > vybound
                    z.zuy[i, j, k] *= vybound / yvel
                end
                if zvel > vzbound
                    z.zuz[i, j, k] *= vzbound / zvel
                end
            end
        end
    end

    return nothing
end
