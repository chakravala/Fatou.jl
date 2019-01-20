#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2019 Michael Reed

using PyPlot
import PyPlot: plot

function plot(K::FilledSet;c::String="",bare::Bool=false)
    # plot figure using imshow based in input preferences
    figure()
    isempty(c) && (c = K.meta.cmap)
    isempty(c) ? PyPlot.imshow(K.meta.iter ? K.iter : K.mix, extent=K.meta.∂) :
        PyPlot.imshow(K.meta.iter ? K.iter : K.mix, cmap=c, extent=K.meta.∂)
    tight_layout()
    if !bare
        # determine if plot is Iteration, Roots, or Limit
        typeof(K.meta.iter ? K.iter : K.mix) == Matrix{UInt8} ? t = L"iter. " :
            K.meta.m==1 ? t = L"roots" : t = L"limit"
        # annotate title using LaTeX
        ttext = "f:z\\mapsto $(rdpm(Algebra.latex(K.meta.E))),\\,"
        if K.meta.newt
            title(latexstring("$ttext m = $(K.meta.m), ")*t)
            # annotate y-axis with Newton's method
            ylabel(L"Fatou\,set:\,"*L"z\,↦\,z-m\,×\,f(z)\,/\,f\,'(z)")
        else
            title(latexstring(ttext)*t)
        end
        tight_layout()
        colorbar()
    end
end

function orbit(E,f::Function,bi::Matrix{Float64},orb::Int=0,depth::Int=1,incr::Int=384; plt::Function=plot)
    x,N,N2,orbit,bis = real_orb(E,f,bi,orb,depth,incr)
    # prepare for next figure
    figure()
    # plot background lines
    plot(x[:],N[:,1],"k--",x[:],N[:,2])
    # plot orbit cobweb path
    plot(orbit[:,1],orbit[:,2],"r")
    # plot f^2,f^3,f^4,...
    for h ∈ 3:depth+1
        plot(x,N[:,h],lw=1)
    end
    if ~(orb == 0)
        plt(range(bi[1],stop=bi[2],length=length(N2)),N2[:],"gray",marker="x",linestyle=":",lw=1)
        funs = [latexstring("\\phi(x_{0:$orb})")]
        funt = ", IC: \$ x_0 = $(bis[3])\$, \$ n\\in0:$orb\$"
    else
        funs = []
        funt = ""
    end
    # trim graph
    d=1.07
    xlim(bi[1],bi[2])
    ylim(minimum([d*minimum(N[:,2]),0]),maximum([d*maximum(N[:,2]),0]))
    # set title
    fune = rdpm(Algebra.latex(E))
    PyPlot.title(latexstring("\$ x \\mapsto $fune\$$funt"))
    # set legend
    legend(vcat([L"$y=x$",L"$\phi(x)$",L"(x_n,\phi(x_n))"],[latexstring("\\phi^{$x}(x)") for x ∈ 2:depth],funs))
    tight_layout()
    return
end
