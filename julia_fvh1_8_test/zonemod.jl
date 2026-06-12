# Julia port of zonemod.f90
# Mirrors the 3D grid and boundary-condition state used by the Fortran code.

module Zone

imax::Int = 40
jmax::Int = 40
kmax::Int = 40
const ndim::Int = 3

zro::Array{Float64,3}  = zeros(Float64, imax, jmax, kmax)
zpr::Array{Float64,3}  = zeros(Float64, imax, jmax, kmax)
zux::Array{Float64,3}  = zeros(Float64, imax, jmax, kmax)
zuy::Array{Float64,3}  = zeros(Float64, imax, jmax, kmax)
zuz::Array{Float64,3}  = zeros(Float64, imax, jmax, kmax)

dpi::Array{Float64,5} = zeros(Float64, 3, 3, imax, jmax, kmax)

repixky::Vector{Float64} = zeros(Float64, jmax)
impixky::Vector{Float64} = zeros(Float64, jmax)
spec::Vector{Float64}    = zeros(Float64, jmax)
spec2::Vector{Float64}   = zeros(Float64, jmax)

rhokx::Vector{Float64} = zeros(Float64, imax)

zxa::Vector{Float64} = zeros(Float64, imax)
zdx::Vector{Float64} = zeros(Float64, imax)
zxc::Vector{Float64} = zeros(Float64, imax)

zya::Vector{Float64} = zeros(Float64, jmax)
zdy::Vector{Float64} = zeros(Float64, jmax)
zyc::Vector{Float64} = zeros(Float64, jmax)

zza::Vector{Float64} = zeros(Float64, kmax)
zdz::Vector{Float64} = zeros(Float64, kmax)
zzc::Vector{Float64} = zeros(Float64, kmax)

nleftx::Int = 0
nlefty::Int = 0
nleftz::Int = 0

nrightx::Int = 0
nrighty::Int = 0
nrightz::Int = 0

end
