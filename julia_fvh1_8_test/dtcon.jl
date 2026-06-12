# Julia port of dtcon.f90
# Limits the velocity components to CFL-style bounds for the fixed timestep.

include("vh1mods.jl")
include("zonemod.jl")

using .Global
using .Zone

"""
    vcon!()

Clamp the velocity components in the Zone arrays to the CFL bounds implied by dt.
"""
function vcon!()
    vxbound = Zone.zdx[1] / Global.dt * 0.45
    vybound = Zone.zdy[1] / Global.dt * 0.45
    vzbound = Zone.zdz[1] / Global.dt * 0.45

    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                xvel = abs(Zone.zux[i, j, k])
                yvel = abs(Zone.zuy[i, j, k])
                zvel = abs(Zone.zuz[i, j, k])

                if xvel > vxbound
                    Zone.zux[i, j, k] *= vxbound / xvel
                end
                if yvel > vybound
                    Zone.zuy[i, j, k] *= vybound / yvel
                end
                if zvel > vzbound
                    Zone.zuz[i, j, k] *= vzbound / zvel
                end
            end
        end
    end

    return nothing
end
