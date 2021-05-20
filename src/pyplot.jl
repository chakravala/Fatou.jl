#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2019 Michael Reed

PyPlot.plot(K::FilledSet;c::String="",bare::Bool=false) = PyPlot.imshow(K;cmap=c,bare=bare)

function PyPlot.imshow(K::FilledSet;cmap::String="",bare::Bool=false)
    PyPlot.figure() # plot figure using imshow based in input preferences
    isempty(cmap) && (cmap = K.meta.cmap)
    isempty(cmap) ? PyPlot.imshow(K.meta.iter ? K.iter : K.mix, extent=K.set.∂) :
        PyPlot.imshow(K.meta.iter ? K.iter : K.mix, cmap=cmap, extent=K.set.∂)
    PyPlot.tight_layout()
    !bare && PyPlot.title(K)
end

function PyPlot.title(K::FilledSet) # annotate title using LaTeX
    text,t = "f:z\\mapsto $(rdpm(Algebra.latex(K.meta.E))),\\,",LaTeXString(typeplot(K))
    if K.meta.newt # annotate y-axis with Newton's method
        PyPlot.title(latexstring("$text m = $(K.meta.m), ")*t)
        PyPlot.ylabel(L"Fatou\,set:\,"*L"z\,↦\,z-m\,×\,f(z)\,/\,f\,'(z)")
    else
        PyPlot.title(latexstring(ttext)*t)
    end
    PyPlot.tight_layout()
    PyPlot.colorbar()
end

function orbit(E,f::Function,bi::Matrix{Float64},orb::Int=0,depth::Int=1,incr::Int=384)
    x,N,N2,orbit,bis = real_orb(E,f,bi,orb,depth,incr)
    # prepare for next figure
    PyPlot.figure()
    # plot background lines
    PyPlot.plot(x[:],N[:,1],"k--",x[:],N[:,2])
    # plot orbit cobweb path
    PyPlot.plot(orbit[:,1],orbit[:,2],"r")
    # plot f^2,f^3,f^4,...
    for h ∈ 3:depth+1
        PyPlot.plot(x,N[:,h],lw=1)
    end
    if ~(orb == 0)
        PyPlot.plot(range(bi[1],stop=bi[2],length=length(N2)),N2[:],"gray",marker="x",linestyle=":",lw=1)
        funs = [latexstring("\\phi(x_{0:$orb})")]
        funt = ", IC: \$ x_0 = $(bis[3])\$, \$ n\\in0:$orb\$"
    else
        funs = []
        funt = ""
    end
    # trim graph
    d=1.07
    PyPlot.xlim(bi[1],bi[2])
    PyPlot.ylim(minimum([d*minimum(N[:,2]),0]),maximum([d*maximum(N[:,2]),0]))
    # set title
    fune = rdpm(Algebra.latex(E))
    PyPlot.title(latexstring("\$ x \\mapsto $fune\$$funt"))
    # set legend
    PyPlot.legend(vcat([L"$y=x$",L"$\phi(x)$",L"(x_n,\phi(x_n))"],[latexstring("\\phi^{$x}(x)") for x ∈ 2:depth],funs))
    PyPlot.tight_layout()
    return
end
