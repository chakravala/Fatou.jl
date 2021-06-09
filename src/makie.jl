#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2019 Michael Reed

import ColorSchemes

AbstractPlotting.plot(K::FilledSet;bare=false,args...) = AbstractPlotting.heatmap(K;bare=bare,args...)

for plt ∈ (:contour,:contourf,:heatmap)
    plt! = Symbol(plt,:!)
    @eval function AbstractPlotting.$plt(K::FilledSet;bare::Bool=false,args...)
        scene, layout = AbstractPlotting.layoutscene()
        ax = layout[1,1] = if K.meta.newt 
            AbstractPlotting.Axis(scene,title=String(K),ylabel="Fatou set: z ↦ z-m×f(z)/f'(z)")
        else
            AbstractPlotting.Axis(scene,title=String(K))
        end
        plt = AbstractPlotting.$plt!(ax,K;args...)
        layout[1,2] = AbstractPlotting.Colorbar(scene,plt;width=30)
        return scene
    end
    @eval function AbstractPlotting.$plt!(ax,K::FilledSet;bare::Bool=false,args...)
        !haskey(args,:colormap) && !isempty(K.meta.cmap) && (return AbstractPlotting.$plt!(ax,K;bare=bare,colormap=Symbol(K.meta.cmap),args...))
        r1,r2 = ranges(K)
        Z = $(plt≠:heatmap ? :(transpose(K.meta.iter ? K.iter : K.mix)) : :(reverse(transpose(K.meta.iter ? K.iter : K.mix),dims=2)))
        if bare
            AbstractPlotting.$plt!(ax,r1,r2,Z;args...)
        else
            AbstractPlotting.$plt!(ax,r1,r2,Z;axis=(title=String(K),),args...)
        end
    end
end # center!(scene)

function AbstractPlotting.surface(K::FilledSet;bare::Bool=false,args...)
    !haskey(args,:colormap) && !isempty(K.meta.cmap) && (return AbstractPlotting.surface(K;bare=bare,colormap=Symbol(K.meta.cmap),args...))
    r1,r2 = ranges(K)
    Z = reverse(reverse(transpose(K.meta.iter ? K.iter : K.mix),dims=1),dims=2)
    if bare
        AbstractPlotting.surface(reverse(r1),r2,Z;(xreversed=true,)args...)
    else
        AbstractPlotting.surface(reverse(r1),r2,Z;axis=(title=String(K),xreversed=true),args...)
    end
end

function AbstractPlotting.arrows(K::FilledSet;bare::Bool=false,args...)
    r1,r2 = ranges(K)
    Z = transpose(K.set.Ω)
    AbstractPlotting.arrows(r1,r2,real.(Z),imag.(Z);args...)
end

function orbit(E,f::Function,bi::Matrix{Float64},orb::Int=0,depth::Int=1,incr::Int=384)
    x,N,N2,orbit,bis = real_orb(E,f,bi,orb,depth,incr)
    # setup plot
    funs,funt = ~(orb == 0) ? (["ϕ(x₀₋ₙ)"],", IC: x₀ = $(bis[3]), n∈0:$orb") : ([],"")
    scene, layout = AbstractPlotting.layoutscene()
    ax = layout[1,1] = AbstractPlotting.Axis(scene,title="x ↦ $E$funt")
    # plot background lines
    plt = [AbstractPlotting.lines!(ax,x[:],N[:,1];linestyle=:dash),
    AbstractPlotting.lines!(ax,x[:],N[:,2]),
    # plot orbit cobweb path
    AbstractPlotting.lines!(ax,orbit[:,1],orbit[:,2];color=:red)]
    # plot f^2,f^3,f^4,...
    for h ∈ 3:depth+1
        push!(plt,AbstractPlotting.lines!(ax,x,N[:,h];linewidth=1,color=ColorSchemes.tab10[h>5 ? h-1 : h-2]))
    end
    if orb≠0
        ran = range(bi[1],stop=bi[2],length=length(N2))
        plt = [plt...,[
            AbstractPlotting.lines!(ax,ran,N2[:];color=:gray,linestyle=:dot,linewidth=1),
            AbstractPlotting.scatter!(ax,ran,N2[:];color=:gray,marker=:x)]]
    end
    # trim graph
    d=1.07
    AbstractPlotting.xlims!(ax,bi[1],bi[2])
    AbstractPlotting.ylims!(ax,minimum([d*minimum(N[:,2]),0]),maximum([d*maximum(N[:,2]),0]))
    # set legend
    layout[1,2] = Legend(scene,plt,vcat(["y=x","ϕ(x)","(xₙ,ϕ(xₙ))"],["ϕ^$x(x)" for x ∈ 2:depth],funs))
    return scene
end

