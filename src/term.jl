#   This file is part of Fatou.jl. It is licensed under the MIT license       
#   Copyright (C) 2019 Michael Reed

using ImageInTerminal, ColorSchemes

nonan(x) = isnan(x) ? 0.0 : x

function Base.show(io::IO,K::FilledSet;c::String="",bare::Bool=false)
    # plot figure using imshow based in input preferences
    isempty(c) && (c = K.meta.cmap)
    S = size(K.iter)
    H = zeros(ColorSchemes.RGB{Float64},S...)
    C = getproperty(ColorSchemes, isempty(c) ? :balance : Symbol(c))
    if K.meta.iter
        M = length(C)/(max(K.iter...)+2)
        for x ∈ 1:S[1], y ∈ 1:S[2]
            H[x,y] = C[round(Int,M*K.iter[x,y]+1)]
        end
    else
        for x ∈ 1:S[1], y ∈ 1:S[2]
            H[x,y] = get(C,nonan(K.mix[x,y]))
        end
    end
    display(H)
    if !bare
        # determine if plot is Iteration, Roots, or Limit
        typeof(K.meta.iter ? K.iter : K.mix) == Matrix{UInt8} ? t = "iter. " :
            K.meta.m==1 ? t = "roots" : t = "limit"
        # annotate title using LaTeX
        ttext = "f : z ↦ $(K.meta.E), "
        if K.meta.newt
            print("$ttext m = $(K.meta.m), ",t)
        else
            print(ttext,t)
        end
    end
end

