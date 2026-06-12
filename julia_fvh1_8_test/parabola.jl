# Julia port of parabola.f90
# Uses the already-ported sweeps-size module from vh1mods.jl.

include("vh1mods.jl")

using .Sweepsize

function parabola!(nmin0::Integer, nmax0::Integer, para::Vector{Float64},
                   a::Vector{Float64}, deltaa::Vector{Float64},
                   a6::Vector{Float64}, al::Vector{Float64},
                   flat::Vector{Float64})
		   
    diffa = zeros(Float64, Sweepsize.maxsweep)
    da = zeros(Float64, Sweepsize.maxsweep)
    ar = zeros(Float64, Sweepsize.maxsweep)
    scrch1 = zeros(Float64, Sweepsize.maxsweep)
    scrch2 = zeros(Float64, Sweepsize.maxsweep)
    scrch3 = zeros(Float64, Sweepsize.maxsweep)

    for n in nmin0:(nmax0 + 1)
        diffa[n] = a[n + 1] - a[n]
    end

    for n in nmin0:nmax0
        da[n] = para[4] * diffa[n] + para[5] * diffa[n - 1]
        da[n] = copysign(min(abs(da[n]), 2.0 * abs(diffa[n - 1]), 2.0 * abs(diffa[n])), da[n])
    end

    for n in nmin0:nmax0
        if diffa[n - 1] * diffa[n] < 0.0
            da[n] = 0.0
        end
    end

    for n in nmin0:nmax0
        ar[n] = a[n] + para[1] * diffa[n] + para[2] * da[n + 1] + para[3] * da[n]
        al[n + 1] = ar[n]
    end

    for n in nmin0:nmax0
        onemfl = 1.0 - flat[n]
        ar[n] = flat[n] * a[n] + onemfl * ar[n]
        al[n] = flat[n] * a[n] + onemfl * al[n]
    end

    for n in nmin0:nmax0
        deltaa[n] = ar[n] - al[n]
        a6[n] = 6.0 * (a[n] - 0.5 * (al[n] + ar[n]))
        scrch1[n] = (ar[n] - a[n]) * (a[n] - al[n])
        scrch2[n] = deltaa[n] * deltaa[n]
        scrch3[n] = deltaa[n] * a6[n]
    end

    for n in nmin0:nmax0
        if scrch1[n] <= 0.0
            ar[n] = a[n]
            al[n] = a[n]
        end
        if scrch2[n] < scrch3[n]
            al[n] = 3.0 * a[n] - 2.0 * ar[n]
        end
        if scrch2[n] < -scrch3[n]
            ar[n] = 3.0 * a[n] - 2.0 * al[n]
        end
    end

    for n in nmin0:nmax0
        deltaa[n] = ar[n] - al[n]
        a6[n] = 6.0 * (a[n] - 0.5 * (al[n] + ar[n]))
    end

    return nothing
end

function paraset!(para::Vector{Float64})
    para[1] = 0.5
    para[2] = -1.0 / 6.0
    para[3] = 1.0 / 6.0
    para[4] = 0.5
    para[5] = 0.5
    return nothing
end
