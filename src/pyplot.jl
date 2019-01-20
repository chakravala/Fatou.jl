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

include("orbitplot.jl")
