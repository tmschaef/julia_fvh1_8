# Constants and state types for the VH-1 Julia port.

using Random

const maxsweep = 1036

const imax = 40
const jmax = 40
const kmax = 40
const ndim = 3

const smallp = 1.0e-15
const smallr = 1.0e-15
const small = 1.0e-15

const pi = 2.0 * asin(1.0)
const gam = 5.0 / 3.0
const gamm = gam - 1.0
const alpha = 0.0
const alphat = 0.0
const nav = 400.0
const tav = 1.0
const etaav = 0.1 * nav

const A0 = 0.25 * 0.25
const V0 = A0 * 0.25

const uinflo = 0.0
const dinflo = 0.0
const vinflo = 0.0
const winflo = 0.0
const pinflo = 0.0
const einflo = 0.0

const uotflo = 0.0
const dotflo = 0.0
const votflo = 0.0
const wotflo = 0.0
const potflo = 0.0
const eotflo = 0.0

mutable struct globalvars
    ncycle::Int
    ncycp::Int
    ntrial::Int
    nacc::Int
    time::Float64
    dt::Float64
    timep::Float64
    dtfix::Float64
    endtime::Float64
    tprin::Float64
end

function globalvars()
    globalvars(0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

mutable struct Sweeps
    nmin::Int
    nmax::Int
    nleft::Int
    nright::Int
    r::Vector{Float64}
    p::Vector{Float64}
    e::Vector{Float64}
    q::Vector{Float64}
    u::Vector{Float64}
    v::Vector{Float64}
    w::Vector{Float64}
    xa::Vector{Float64}
    xa0::Vector{Float64}
    dx::Vector{Float64}
    dx0::Vector{Float64}
    dvol::Vector{Float64}
    dvol0::Vector{Float64}
    flat::Vector{Float64}
    radius::Float64
    area::Float64
end

function Sweeps()
    Sweeps(
        0, 0, 0, 0,
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        zeros(Float64, maxsweep),
        0.0, 0.0
    )
end

mutable struct zone
    zro::Array{Float64,3}
    zpr::Array{Float64,3}
    zux::Array{Float64,3}
    zuy::Array{Float64,3}
    zuz::Array{Float64,3}
end

function zone()
    zone(
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax)
    )
end

mutable struct zonegrid
    zxa::Vector{Float64}
    zdx::Vector{Float64}
    zxc::Vector{Float64}
    zya::Vector{Float64}
    zdy::Vector{Float64}
    zyc::Vector{Float64}
    zza::Vector{Float64}
    zdz::Vector{Float64}
    zzc::Vector{Float64}
    nleftx::Int
    nlefty::Int
    nleftz::Int
    nrightx::Int
    nrighty::Int
    nrightz::Int
end

function zonegrid()
    zonegrid(
        zeros(Float64, imax),
        zeros(Float64, imax),
        zeros(Float64, imax),
        zeros(Float64, jmax),
        zeros(Float64, jmax),
        zeros(Float64, jmax),
        zeros(Float64, kmax),
        zeros(Float64, kmax),
        zeros(Float64, kmax),
        0, 0, 0, 0, 0, 0
    )
end

mutable struct zoneobs
    repixky::Vector{Float64}
    impixky::Vector{Float64}
    spec::Vector{Float64}
    spec2::Vector{Float64}
    rhokx::Vector{Float64}
end

function zoneobs()
    zoneobs(
        zeros(Float64, jmax),
        zeros(Float64, jmax),
        zeros(Float64, jmax),
        zeros(Float64, jmax),
        zeros(Float64, imax)
    )
end

# Reusable scratch buffers for ppmlr / parabola / remap / evolve / states.
mutable struct PpmWork
    para::Vector{Float64}
    diffa::Vector{Float64}
    da::Vector{Float64}
    ar::Vector{Float64}
    scrch1::Vector{Float64}
    scrch2::Vector{Float64}
    scrch3::Vector{Float64}
    dr::Vector{Float64}
    du::Vector{Float64}
    dp::Vector{Float64}
    r6::Vector{Float64}
    u6::Vector{Float64}
    p6::Vector{Float64}
    rl::Vector{Float64}
    ul::Vector{Float64}
    pl::Vector{Float64}
    rrgh::Vector{Float64}
    urgh::Vector{Float64}
    prgh::Vector{Float64}
    rlft::Vector{Float64}
    ulft::Vector{Float64}
    plft::Vector{Float64}
    umid::Vector{Float64}
    pmid::Vector{Float64}
    Cdtdx::Vector{Float64}
    fCdtdx::Vector{Float64}
    amid::Vector{Float64}
    uold::Vector{Float64}
    xa1::Vector{Float64}
    dvol1::Vector{Float64}
    upmid::Vector{Float64}
    dm::Vector{Float64}
    dtbdm::Vector{Float64}
    xa2::Vector{Float64}
    xa3::Vector{Float64}
    dv::Vector{Float64}
    vl::Vector{Float64}
    v6::Vector{Float64}
    dw::Vector{Float64}
    wl::Vector{Float64}
    w6::Vector{Float64}
    de::Vector{Float64}
    el::Vector{Float64}
    e6::Vector{Float64}
    dq::Vector{Float64}
    ql::Vector{Float64}
    q6::Vector{Float64}
    dm0::Vector{Float64}
    delta::Vector{Float64}
    fluxr::Vector{Float64}
    fluxu::Vector{Float64}
    fluxv::Vector{Float64}
    fluxw::Vector{Float64}
    fluxe::Vector{Float64}
    fluxq::Vector{Float64}
end

function _ppm_vec()
    zeros(Float64, maxsweep)
end

function PpmWork()
    PpmWork(
        zeros(Float64, 5),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(),
        _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec(), _ppm_vec()
    )
end

# Shared grid buffers for mc3d (read/written by non-overlapping cell sets per color phase).
mutable struct McWork
    p1::Array{Float64,3}
    p2::Array{Float64,3}
    p3::Array{Float64,3}
    q::Array{Float64,3}
end

function McWork()
    McWork(
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax),
        zeros(Float64, imax, jmax, kmax)
    )
end

# Per-thread scratch + RNG for parallel mc3d (no per-cell allocation).
mutable struct McScratch
    lambda::Matrix{Float64}
    lambda_t::Matrix{Float64}
    delta_p::Matrix{Float64}
    qnew_cube::Vector{Float64}
    pcube1::Vector{Float64}
    pcube2::Vector{Float64}
    pcube3::Vector{Float64}
    rng::Xoshiro
end

function McScratch(tid::Int)
    McScratch(
        zeros(Float64, 3, 3),
        zeros(Float64, 3, 3),
        zeros(Float64, 3, 3),
        zeros(Float64, 8),
        zeros(Float64, 8),
        zeros(Float64, 8),
        zeros(Float64, 8),
        Xoshiro(123456 + tid * 100003)
    )
end

mutable struct McThreadPool
    grid::McWork
    scratch::Vector{McScratch}
end

function McThreadPool()
    n = Threads.maxthreadid()
    McThreadPool(McWork(), [McScratch(tid) for tid in 1:n])
end

@inline mc_phase_count(n::Int, parity::Int) = (n - parity + 1) ÷ 2

# Per-thread Sweeps + PpmWork for parallel line sweeps (sweepx/y/z).
mutable struct SweepThreadPool
    sweeps::Vector{Sweeps}
    work::Vector{PpmWork}
end

function SweepThreadPool()
    n = Threads.maxthreadid()
    SweepThreadPool([Sweeps() for _ in 1:n], [PpmWork() for _ in 1:n])
end

function _build_pixky_trig()
    nk = jmax
    cos_j = Matrix{Float64}(undef, jmax, nk)
    sin_j = Matrix{Float64}(undef, jmax, nk)
    sin_i = Matrix{Float64}(undef, imax, nk)
    for in_ in 1:nk
        for j in 1:jmax
            ang = 2.0 * π * j * in_ / jmax
            cos_j[j, in_] = cos(ang)
            sin_j[j, in_] = sin(ang)
        end
        for i in 1:imax
            sin_i[i, in_] = sin(2.0 * π * (i - 1) * in_ / imax)
        end
    end
    return cos_j, sin_j, sin_i
end

const PIXKY_COS_J, PIXKY_SIN_J, PIXKY_SIN_I = _build_pixky_trig()
const PIXKY_XNORM = 1.0 / sqrt(Float64(imax * jmax * kmax))

const MC_SIGN_I = (-1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0)
const MC_SIGN_J = (-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0)
const MC_SIGN_K = (-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0)
const MC_I_SLOT = (1, 2, 1, 2, 1, 2, 1, 2)
const MC_J_SLOT = (1, 1, 2, 2, 1, 1, 2, 2)
const MC_K_SLOT = (1, 1, 1, 1, 2, 2, 2, 2)
