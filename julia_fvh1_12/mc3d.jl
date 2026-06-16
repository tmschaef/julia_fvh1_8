# Julia port of mc3d.f90 — 8-color checkerboard parallelization.

@inline function gaussian_random(rng::Random.AbstractRNG)
    r1 = rand(rng)
    r2 = rand(rng)
    if r1 <= 0.0
        r1 = 1.0e-18
    end
    return sqrt(-2.0 * log(r1)) * cos(2.0 * π * r2)
end

@inline function mc_corner_update!(dhtot::Float64, qnew_cube, pcube1, pcube2, pcube3,
                                   p1, p2, p3, q, zro,
                                   ci::Int, cj::Int, ck::Int,
                                   si::Float64, sj::Float64, sk::Float64,
                                   dp11, dp12, dp13, dp21, dp22, dp23, dp31, dp32, dp33,
                                   cube_idx::Int)
    p_old1 = p1[ci, cj, ck]
    p_old2 = p2[ci, cj, ck]
    p_old3 = p3[ci, cj, ck]

    p_new1 = p_old1 + 0.25 * (si * dp11 + sj * dp12 + sk * dp13)
    p_new2 = p_old2 + 0.25 * (si * dp21 + sj * dp22 + sk * dp23)
    p_new3 = p_old3 + 0.25 * (si * dp31 + sj * dp32 + sk * dp33)

    pcube1[cube_idx] = p_new1
    pcube2[cube_idx] = p_new2
    pcube3[cube_idx] = p_new3

    mass = V0 * zro[ci, cj, ck]
    H_old = (p_old1^2 + p_old2^2 + p_old3^2) / (2.0 * mass)
    H_new = (p_new1^2 + p_new2^2 + p_new3^2) / (2.0 * mass)

    dhtot += H_new - H_old
    qnew_cube[cube_idx] = q[ci, cj, ck] - (H_new - H_old) / mass
    return dhtot
end

@inline function mc3d_cell!(rng::Random.AbstractRNG,
                            p1, p2, p3, q, zro,
                            sc::McScratch,
                            i::Int, j::Int, k::Int,
                            factor::Float64)
    lambda = sc.lambda
    lambda_t = sc.lambda_t
    delta_p = sc.delta_p
    qnew_cube = sc.qnew_cube
    pcube1 = sc.pcube1
    pcube2 = sc.pcube2
    pcube3 = sc.pcube3

    i_corners = (mod(i - 1, imax) + 1, mod(i, imax) + 1)
    j_corners = (mod(j - 1, jmax) + 1, mod(j, jmax) + 1)
    k_corners = (mod(k - 1, kmax) + 1, mod(k, kmax) + 1)

    @inbounds for col in 1:3, row in 1:3
        lambda[row, col] = gaussian_random(rng)
    end

    trace_lambda = lambda[1, 1] + lambda[2, 2] + lambda[3, 3]
    @inbounds for col in 1:3, row in 1:3
        lambda_t[row, col] = lambda[row, col] + lambda[col, row]
    end
    lambda_t[1, 1] -= 2.0 / 3.0 * trace_lambda
    lambda_t[2, 2] -= 2.0 / 3.0 * trace_lambda
    lambda_t[3, 3] -= 2.0 / 3.0 * trace_lambda
    @inbounds for col in 1:3, row in 1:3
        delta_p[row, col] = factor * lambda_t[row, col]
    end

    dp11 = delta_p[1, 1]
    dp12 = delta_p[1, 2]
    dp13 = delta_p[1, 3]
    dp21 = delta_p[2, 1]
    dp22 = delta_p[2, 2]
    dp23 = delta_p[2, 3]
    dp31 = delta_p[3, 1]
    dp32 = delta_p[3, 2]
    dp33 = delta_p[3, 3]

    dhtot = 0.0
    @inbounds for cube_idx in 1:8
        ci = i_corners[MC_I_SLOT[cube_idx]]
        cj = j_corners[MC_J_SLOT[cube_idx]]
        ck = k_corners[MC_K_SLOT[cube_idx]]
        dhtot = mc_corner_update!(
            dhtot, qnew_cube, pcube1, pcube2, pcube3,
            p1, p2, p3, q, zro,
            ci, cj, ck,
            MC_SIGN_I[cube_idx], MC_SIGN_J[cube_idx], MC_SIGN_K[cube_idx],
            dp11, dp12, dp13, dp21, dp22, dp23, dp31, dp32, dp33,
            cube_idx
        )
    end

    accept = false
    if rand(rng) < exp(-dhtot / tav)
        accept = true
        @inbounds for cube_idx in 1:8
            if qnew_cube[cube_idx] <= 0.0
                accept = false
                break
            end
        end
    end

    if accept
        @inbounds for cube_idx in 1:8
            ci = i_corners[MC_I_SLOT[cube_idx]]
            cj = j_corners[MC_J_SLOT[cube_idx]]
            ck = k_corners[MC_K_SLOT[cube_idx]]
            p1[ci, cj, ck] = pcube1[cube_idx]
            p2[ci, cj, ck] = pcube2[cube_idx]
            p3[ci, cj, ck] = pcube3[cube_idx]
            q[ci, cj, ck] = qnew_cube[cube_idx]
        end
        return true
    end
    return false
end

function mc3d!(g::globalvars, z::zone, pool::McThreadPool)
    mc = pool.grid
    p1 = mc.p1
    p2 = mc.p2
    p3 = mc.p3
    q = mc.q

    p1 .= V0 .* z.zro .* z.zux
    p2 .= V0 .* z.zro .* z.zuy
    p3 .= V0 .* z.zro .* z.zuz
    q  .= z.zpr ./ (gamm .* z.zro)

    factor = A0 * sqrt(etaav * tav * g.dt / V0)

    ntrial_add = zeros(Int, Threads.maxthreadid())
    nacc_add = zeros(Int, Threads.maxthreadid())

    # 8 phases: cells with (i%2, j%2, k%2) == (pi, pj, pk).
    # Same-phase cells are at least Chebyshev distance 2 apart, so corner
    # writes do not overlap within a phase.
    @inbounds for phase in 0:7
        pi = phase & 1
        pj = (phase >> 1) & 1
        pk = (phase >> 2) & 1

        ni = mc_phase_count(imax, pi)
        nj = mc_phase_count(jmax, pj)
        nk = mc_phase_count(kmax, pk)
        ncells = ni * nj * nk

        Threads.@threads for idx in 1:ncells
            tid = Threads.threadid()
            sc = pool.scratch[tid]
            rng = sc.rng

            tmp = idx - 1
            ik = tmp % nk + 1
            tmp = tmp ÷ nk
            ij = tmp % nj + 1
            ii = tmp ÷ nj + 1

            i = 2 * (ii - 1) + (pi == 0 ? 2 : 1)
            j = 2 * (ij - 1) + (pj == 0 ? 2 : 1)
            k = 2 * (ik - 1) + (pk == 0 ? 2 : 1)

            ntrial_add[tid] += 1
            if mc3d_cell!(rng, p1, p2, p3, q, z.zro, sc, i, j, k, factor)
                nacc_add[tid] += 1
            end
        end
    end

    g.ntrial += sum(ntrial_add)
    g.nacc += sum(nacc_add)

    inv_v0 = 1.0 / V0
    @inbounds for k in 1:kmax, j in 1:jmax, i in 1:imax
        rho_inv = inv_v0 / z.zro[i, j, k]
        z.zux[i, j, k] = p1[i, j, k] * rho_inv
        z.zuy[i, j, k] = p2[i, j, k] * rho_inv
        z.zuz[i, j, k] = p3[i, j, k] * rho_inv
        z.zpr[i, j, k] = gamm * q[i, j, k] * z.zro[i, j, k]
    end

    return nothing
end
