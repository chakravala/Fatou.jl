#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export orbit

"""
    orbit(K::Fatou.Define)
    orbit(u::Function, ∂, orbit, depth, n)

Plot funciton compositions of the primary `Fatou.Define` function up to any `depth` including cobweb plot with an `orbit` depth from the `x0` start point.

# Examples
```Julia
julia> juliafill("z^2-0.67",∂=[-1.25,1.5],x0=1.25,orbit=17,depth=3) |> orbit
```
"""
function orbit(K::Define{FT,QT,CT}) where {FT,QT,CT}
    !isdefined(Fatou,:UnicodePlots) && !isdefined(Fatou,:PyPlot) && throw(error("Requires `using PyPlot` or `using UnicodePlots`"))
    K.x0 == nothing ? (bi = K.∂[1:2]') : (bi = [K.∂[1:2]...,K.x0]')
    orbit(K.E,z->K.F(z,0),convert(Array{Float64},bi),K.orbit,K.depth,Int(K.n))
end

function real_orb(E,f::Function,bi::Matrix{Float64},orb::Int=0,depth::Int=1,incr::Int=384)
    # initalize array to depth
    N = zeros(incr,depth+1)
    bis = zeros(3)
    bis[1:length(bi)] = bi[:]
    # set x-axis coordinate set
    x = range(bi[1],stop=bi[2],length=incr)
    N[:,1] = x[:]
    # loop over all discrete x-axis points in set
    # loop function composition at x val
    for t ∈ 1:depth
        N[:,t+1] = broadcast(f,N[:,t])
    end
    # initialize for orbit cobweb
    N2 = zeros(orb+1)
    N2[1] = bis[3];
    # loop over orbit cobweb
    for t ∈ 1:orb
        # evaluate function composition
        N2[t+1] = f(N2[t])
    end
    # interleave cobweb orbit data
    siz = 3*(length(N2)-1)
    orbit = zeros(siz,2)
    orbit[1:3:siz,1] = N2[1:end-1]
    orbit[2:3:siz,1] = N2[1:end-1]
    orbit[3:3:siz,1] = N2[2:end]
    orbit[1:3:siz,2] = N2[1:end-1]
    orbit[2:3:siz,2] = N2[2:end]
    orbit[3:3:siz,2] = N2[2:end]
    return x,N,N2,orbit,bis
end
