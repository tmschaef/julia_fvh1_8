# Julia port of mc3d.f90
# Depends on the already-ported modules in vh1mods.jl and zonemod.jl.

include("vh1mods.jl")
include("zonemod.jl")

using .Global
using .Zone: zro, zpr, zux, zuy, zuz

function wrap_index(i::Integer, n::Integer)
    return mod(i - 1, n) + 1
end

function gaussian_random()
    r1 = rand()
    r2 = rand()
    if r1 <= 0.0
        r1 = 1.0e-18
    end
    return sqrt(-2.0 * log(r1)) * cos(2.0 * π * r2)
end

function mc3d!()
    # Local work arrays
    p1 = zeros(Float64, Zone.imax, Zone.jmax, Zone.kmax)
    p2 = zeros(Float64, Zone.imax, Zone.jmax, Zone.kmax)
    p3 = zeros(Float64, Zone.imax, Zone.jmax, Zone.kmax)
    q  = zeros(Float64, Zone.imax, Zone.jmax, Zone.kmax)

    p1 .= Global.V0 .* Zone.zro .* Zone.zux
    p2 .= Global.V0 .* Zone.zro .* Zone.zuy
    p3 .= Global.V0 .* Zone.zro .* Zone.zuz
    q  .= Zone.zpr ./ (Global.gamm .* Zone.zro)

    factor = Global.A0 * sqrt(Global.etaav * Global.tav * Global.dt / Global.V0)

    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                iidx = zeros(Int, 2)
                jidx = zeros(Int, 2)
                kidx = zeros(Int, 2)

                for ii_idx in 1:2
                    iidx[ii_idx] = wrap_index(i + ii_idx - 1, Zone.imax)
                end
                for jj_idx in 1:2
                    jidx[jj_idx] = wrap_index(j + jj_idx - 1, Zone.jmax)
                end
                for kk_idx in 1:2
                    kidx[kk_idx] = wrap_index(k + kk_idx - 1, Zone.kmax)
                end

                lambda = zeros(Float64, 3, 3)
                for ii_idx in 1:3
                    for jj_idx in 1:3
                        lambda[ii_idx, jj_idx] = gaussian_random()
                    end
                end

                trace_lambda = lambda[1, 1] + lambda[2, 2] + lambda[3, 3]
                lambda_t = lambda + transpose(lambda)
                lambda_t[1, 1] -= 2.0 / 3.0 * trace_lambda
                lambda_t[2, 2] -= 2.0 / 3.0 * trace_lambda
                lambda_t[3, 3] -= 2.0 / 3.0 * trace_lambda
                delta_p = factor .* lambda_t

                cube_i = zeros(Int, 8)
                cube_j = zeros(Int, 8)
                cube_k = zeros(Int, 8)
                sign_i = zeros(Float64, 8)
                sign_j = zeros(Float64, 8)
                sign_k = zeros(Float64, 8)

                cube_idx = 0
                for kk_idx in 1:2
                    for jj_idx in 1:2
                        for ii_idx in 1:2
                            cube_idx += 1
                            cube_i[cube_idx] = iidx[ii_idx]
                            cube_j[cube_idx] = jidx[jj_idx]
                            cube_k[cube_idx] = kidx[kk_idx]
                            sign_i[cube_idx] = (ii_idx == 1) ? -1.0 : 1.0
                            sign_j[cube_idx] = (jj_idx == 1) ? -1.0 : 1.0
                            sign_k[cube_idx] = (kk_idx == 1) ? -1.0 : 1.0
                        end
                    end
                end

                dhtot = 0.0
                qnew_cube = zeros(Float64, 8)
                pcube1 = zeros(Float64, 8)
                pcube2 = zeros(Float64, 8)
                pcube3 = zeros(Float64, 8)

                for cube_idx in 1:8
                    p_old1 = p1[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]]
                    p_old2 = p2[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]]
                    p_old3 = p3[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]]

                    p_new1 = p_old1 + 0.25 * (sign_i[cube_idx] * delta_p[1, 1] +
                                             sign_j[cube_idx] * delta_p[1, 2] +
                                             sign_k[cube_idx] * delta_p[1, 3])
                    p_new2 = p_old2 + 0.25 * (sign_i[cube_idx] * delta_p[2, 1] +
                                             sign_j[cube_idx] * delta_p[2, 2] +
                                             sign_k[cube_idx] * delta_p[2, 3])
                    p_new3 = p_old3 + 0.25 * (sign_i[cube_idx] * delta_p[3, 1] +
                                             sign_j[cube_idx] * delta_p[3, 2] +
                                             sign_k[cube_idx] * delta_p[3, 3])

                    pcube1[cube_idx] = p_new1
                    pcube2[cube_idx] = p_new2
                    pcube3[cube_idx] = p_new3

                    mass = Global.V0 * Zone.zro[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]]
                    H_old = (p_old1^2 + p_old2^2 + p_old3^2) / (2.0 * mass)
                    H_new = (p_new1^2 + p_new2^2 + p_new3^2) / (2.0 * mass)

                    dhtot += H_new - H_old
                    qnew_cube[cube_idx] = q[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]] - (H_new - H_old) / mass
                end

                Global.ntrial += 1
                r = rand()
                accept = false
                if r < exp(-dhtot / Global.tav)
                    accept = true
                    for cube_idx in 1:8
                        if qnew_cube[cube_idx] <= 0.0
                            accept = false
                            break
                        end
                    end
                end

                if accept
                    Global.nacc += 1
                    for cube_idx in 1:8
                        p1[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]] = pcube1[cube_idx]
                        p2[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]] = pcube2[cube_idx]
                        p3[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]] = pcube3[cube_idx]
                        q[cube_i[cube_idx], cube_j[cube_idx], cube_k[cube_idx]] = qnew_cube[cube_idx]
                    end
                end
            end
        end
    end

    for k in 1:Zone.kmax
        for j in 1:Zone.jmax
            for i in 1:Zone.imax
                Zone.zux[i, j, k] = p1[i, j, k] / (Global.V0 * Zone.zro[i, j, k])
                Zone.zuy[i, j, k] = p2[i, j, k] / (Global.V0 * Zone.zro[i, j, k])
                Zone.zuz[i, j, k] = p3[i, j, k] / (Global.V0 * Zone.zro[i, j, k])
                Zone.zpr[i, j, k] = Global.gamm * q[i, j, k] * Zone.zro[i, j, k]
            end
        end
    end

    return nothing
end
