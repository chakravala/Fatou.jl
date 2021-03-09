#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2019 Michael Reed

function orbit(E,f::Function,bi::Matrix{Float64},orb::Int=0,depth::Int=1,incr::Int=384)
    x,N,N2,orbit,bis = real_orb(E,f,bi,orb,depth,incr)
    # trim graph
    d=1.07
    xl = bi[1:2]
    yl = [minimum([d*minimum(N[:,2]),0]),maximum([d*maximum(N[:,2]),0])]
    # plot background lines
    plt = UnicodePlots.lineplot(x[:],N[:,1],xlim=xl,ylim=yl,color=:magenta,name="y=real(z)")
    UnicodePlots.lineplot!(plt,x[:],N[:,2],color=:blue,name="ϕ(real(z))")
    # plot orbit cobweb path
    UnicodePlots.scatterplot!(plt,orbit[:,1],orbit[:,2],color=:red,name="(zₙ,ϕ(zₙ))")
    # set legend
    leg = ["ϕ^{$x}(z)" for x ∈ 2:depth]
    # plot f^2,f^3,f^4,...
    for h ∈ 3:depth+1
        UnicodePlots.lineplot!(plt,x,N[:,h],color=isodd(h) ? :yellow : :green,name=leg[h-2])
    end
    if ~(orb == 0)
        UnicodePlots.lineplot!(plt,collect(range(bi[1],stop=bi[2],length=length(N2))),N2[:],color=:cyan,name="ϕ(z_{0:$orb})")
        funt = ", IC: z₀ = $(bis[3]), n∈0:$orb"
    else
        funt = ""
    end
    # set title
    UnicodePlots.title!(plt,"z ↦ $E$funt")
    return plt
end
