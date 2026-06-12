# Julia port of vh1mods.f90
# Mirrors the global variables and sweep dimensions used by the Fortran code.

module Global

prefix::String = ""                # prefix for output filenames

ncycle::Int = 0
ncycp::Int = 0
ntrial::Int = 0
nacc::Int = 0

time::Float64 = 0.0
 dt::Float64 = 0.0
 timep::Float64 = 0.0
 dtfix::Float64 = 0.0
 endtime::Float64 = 0.0
 tprin::Float64 = 0.0

gam::Float64 = 0.0
pi::Float64 = 0.0
gamm::Float64 = 0.0

alpha::Float64 = 0.0
alphat::Float64 = 0.0

nav::Float64 = 0.0
tav::Float64 = 0.0
etaav::Float64 = 0.0

A0::Float64 = 0.0
V0::Float64 = 0.0

const courant::Float64 = 0.5
const smallp::Float64 = 1.0e-15
const smallr::Float64 = 1.0e-15
const small::Float64  = 1.0e-15

uinflo::Float64 = 0.0
dinflo::Float64 = 0.0
vinflo::Float64 = 0.0
winflo::Float64 = 0.0
pinflo::Float64 = 0.0
einflo::Float64 = 0.0

uotflo::Float64 = 0.0
dotflo::Float64 = 0.0
votflo::Float64 = 0.0
wotflo::Float64 = 0.0
potflo::Float64 = 0.0
eotflo::Float64 = 0.0

end


module Sweepsize

const maxsweep::Int = 1036

end


module Sweeps

using Main.Sweepsize

nmin::Int = 0
nmax::Int = 0
nleft::Int = 0
nright::Int = 0

r::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
p::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
e::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
q::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
u::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
v::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
w::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)

fdiss::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
ediss::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)

xa::Vector{Float64}   = zeros(Float64, Sweepsize.maxsweep)
xa0::Vector{Float64}  = zeros(Float64, Sweepsize.maxsweep)
dx::Vector{Float64}   = zeros(Float64, Sweepsize.maxsweep)
dx0::Vector{Float64}  = zeros(Float64, Sweepsize.maxsweep)
dvol::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)
dvol0::Vector{Float64}= zeros(Float64, Sweepsize.maxsweep)

flat::Vector{Float64} = zeros(Float64, Sweepsize.maxsweep)

radius::Float64 = 0.0
area::Float64 = 0.0

end
