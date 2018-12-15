module Fatou
using SyntaxTree,Reduce,PyPlot,Base.Threads

#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export fatou, juliafill, mandelbrot, newton, basin, plot

"""
    Fatou.Define(E::Any;                  # primary map, (z, c) -> F
      Q::Expr     = :(abs2(z)),           # escape criterion, (z, c) -> Q
      C::Expr     = :((angle(z)/(2π))*n^p)# coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # horizontal grid points
      N::Integer  = 35,                   # max. iterations
      ϵ::Number   = 4,                    # basin ϵ-Limit criterion
      iter::Bool  = false,                # toggle iteration mode
      p::Number   = 0,                    # iteration color exponent
      newt::Bool  = false,                # toggle Newton mode
      m::Number   = 0,                    # Newton multiplicity factor
      O::String   = F,                    # original Newton map
      mandel::Bool= false,                # toggle Mandelbrot mode
      seed::Number= 0.0+0.0im,            # Mandelbrot seed value
      x0          = nothing,              # orbit starting point
      orbit::Int  = 0,                    # orbit cobweb depth
      depth::Int  = 1,                    # depth of function composition
      cmap::String= "")                   # imshow color map

`Define` the metadata for a `Fatou.FilledSet`.
"""
struct Define{FT<:Function,QT<:Function,CT<:Function,M,N}
    E::Any # input expression
    F::FT # primary map
    Q::QT # escape criterion
    C::CT # complex fixed point coloring
    ∂::Array{Float64,1} # bounds
    n::UInt16 # number of grid points
    N::UInt8 # number of iterations
    ϵ::Float64 # epsilon Limit criterion
    iter::Bool # toggle iteration mode
    p::Float64 # iteration color exponent
    newt::Bool # toggle Newton mode
    m::Number # newton multiplicity factor
    mandel::Bool # toggle Mandelbrot mode
    seed::Number # Mandelbrot seed value
    x0 # orbit starting point
    orbit::Int # orbit cobweb depth
    depth::Int # depth of function composition
    cmap::String # imshow color map
    function Define(E;
            Q=:(abs2(z)),
            C=:((angle(z)/(2π))*n^p),
            ∂=π/2,
            n::Integer=176,
            N::Integer=35,
            ϵ::Number=4,
            iter::Bool=false,
            p::Number=0,
            newt::Bool=false,
            m::Number=0,
            mandel::Bool=false,
            seed::Number=0.0+0.0im,
            x0=nothing,
            orbit::Int=0,
            depth::Int=1,
            cmap::String="")
        !(typeof(∂) <: Array) && (∂ = [-float(∂),∂,-∂,∂])
        length(∂) == 2 && (∂ = [∂[1],∂[2],∂[1],∂[2]])
        !newt ? (f = genfun(E,[:z,:c]); q = genfun(Q,[:z,:c])) :
        (f = genfun(newton_raphson(E,m),[:z,:c]); q = genfun(Expr(:call,:abs,E),[:z,:c]))
        c = genfun(C,[:z,:n,:p])
        e = typeof(E) == String ? parse(E) : E
        return new{typeof(f),typeof(q),typeof(c),mandel,newt}(e,f,q,c,convert(Array{Float64,1},∂),UInt16(n),UInt8(N),float(ϵ),iter,float(p),newt,m,mandel,seed,x0,orbit,depth,cmap)
    end
end

"""
    Fatou.FilledSet(::Fatou.Define)

Compute the `Fatou.FilledSet` set using `Fatou.Define`.
"""
struct FilledSet{FT,QT,CT,M,N}
    meta::Define{FT,QT,CT,M,N}
    set::Matrix{Complex{Float64}}
    iter::Matrix{UInt8}
    mix::Matrix{Float64}
    function FilledSet{FT,QT,CT,M,N}(K::Define{FT,QT,CT,M,N}) where {FT,QT,CT,M,N}
        (i,s) = Compute(K)
        return new{FT,QT,CT,M,N}(K,s,i,broadcast(K.C,s,broadcast(float,i./K.N),K.p))
    end
end

  """
      fatou(::Fatou.Define)

Compute the `Fatou.FilledSet` set using `Fatou.Define`.

# Examples
```Julia
julia> fatou(K)
```
"""
fatou(K::Define{FT,QT,CT,M,N}) where {FT,QT,CT,M,N} = FilledSet{FT,QT,CT,M,N}(K)

"""
    juliafill(::Expr;                     # primary map, (z, c) -> F
      Q::Expr     = :(abs2(z)),           # escape criterion, (z, c) -> Q
      C::Expr     = :((angle(z)/(2π))*n^p)# coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # horizontal grid points
      N::Integer  = 35,                   # max. iterations
      ϵ::Number   = 4,                    # basin ϵ-Limit criterion
      iter::Bool  = false,                # toggle iteration mode
      p::Number   = 0,                    # iteration color exponent
      x0          = nothing,              # orbit starting point
      orbit::Int  = 0,                    # orbit cobweb depth
      depth::Int  = 1,                    # depth of function composition
      cmap::String= "")                   # imshow color map

`Define` filled Julia basin in `Fatou`

# Exmaples
```Julia
julia> juliafill("z^2-0.06+0.67im",∂=[-1.5,1.5,-1,1],N=80,n=1501,cmap="RdGy")
```
"""
function juliafill(E;
        Q=:(abs2(z)),
        C=:((angle(z)/(2π))*n^p),
        ∂=π/2,
        n::Integer=176,
        N::Integer=35,
        ϵ::Number=4,
        iter::Bool=false,
        p::Number=0,
        newt::Bool=false,
        m::Number=0,
        x0=nothing,
        orbit::Int=0,
        depth::Int=1,
        cmap::String="")
    return Define(E,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=newt,m=m,x0=x0,orbit=orbit,depth=depth,cmap=cmap)
end

"""
    mandelbrot(::Expr;                    # primary map, (z, c) -> F
      Q::Expr     = :(abs2(z)),           # escape criterion, (z, c) -> Q
      C::Expr     = :(exp(-abs(z))*n^p),  # coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # horizontal grid points
      N::Integer  = 35,                   # max. iterations
      ϵ::Number   = 4,                    # basin ϵ-Limit criterion
      iter::Bool  = false,                # toggle iteration mode
      p::Number   = 0,                    # iteration color exponent
      m::Number   = 0,                    # Newton multiplicity factor
      seed::Number= 0.0+0.0im,            # Mandelbrot seed value
      x0          = nothing,              # orbit starting point
      orbit::Int  = 0,                    # orbit cobweb depth
      depth::Int  = 1,                    # depth of function composition
      cmap::String= "")                   # imshow color map

`Define` Mandelbrot basin in `Fatou`

# Examples
```Julia
mandelbrot(:(z^2+c),n=800,N=20,∂=[-1.91,0.51,-1.21,1.21],cmap="nipy_spectral")
```
"""
function mandelbrot(E;
        Q=:(abs2(z)),
        ∂=π/2,
        C=:(exp(-abs(z))*n^p),
        n::Integer=176,
        N::Integer=35,
        ϵ::Number=4,
        iter::Bool=false,
        p::Number=0,
        newt::Bool=false,
        m::Number=0,
        seed::Number=0.0+0.0im,
        x0=nothing,
        orbit::Int=0,
        depth::Int=1,
        cmap::String="")
    m ≠ 0 && (newt = true)
    return Define(E,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=newt,m=m,mandel=true,seed=seed,x0=x0,orbit=orbit,depth=depth,cmap=cmap)
end

"""
    newton(::Expr;                        # primary map, (z, c) -> F
      C::Expr     = :((angle(z)/(2π))*n^p)# coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # horizontal grid points
      N::Integer  = 35,                   # max. iterations
      ϵ::Number   = 4,                    # basin ϵ-Limit criterion
      iter::Bool  = false,                # toggle iteration mode
      p::Number   = 0,                    # iteration color exponent
      m::Number   = 0,                    # Newton multiplicity factor
      mandel::Bool= false,                # toggle Mandelbrot mode
      seed::Number= 0.0+0.0im,            # Mandelbrot seed value
      x0          = nothing,              # orbit starting point
      orbit::Int  = 0,                    # orbit cobweb depth
      depth::Int  = 1,                    # depth of function composition
      cmap::String= "")                   # imshow color map

`Define` Newton basin in `Fatou`

# Examples
```Julia
julia> newton("z^3-1",n=800,cmap="brg")
```
"""
function newton(E;
        C=:((angle(z)/(2π))*n^p),
        ∂=π/2,
        n::Integer=176,
        N::Integer=35,
        ϵ::Number=0.01,
        iter::Bool=false,
        p::Number=0,
        m::Number=1,
        mandel::Bool=false,
        seed::Number=0.0+0.0im,
        x0=nothing,
        orbit::Int=0,
        depth::Int=1,
        cmap::String="")
    return Define(E,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=true,m=m,mandel=mandel,seed=seed,x0=x0,orbit=orbit,depth=depth,cmap=cmap)
end

# load additional functionality
include("internals.jl"); include("orbitplot.jl")

"""
    basin(::Fatou.Define, ::Integer)

Output the `j`-th basin of `Fatou.Define` as LaTeX.
Each subsequent iteration of the Newton-Raphson method will yield a more complicated set.

# Examples
```Julia
julia> basin(newton("z^3-1"),2)
L"\$\\displaystyle D_2(\\epsilon) = \\left\\{z\\in\\mathbb{C}:\\left|\\,z - \\frac{\\left(z - \\frac{z^{3} - 1}{3 z^{2}}\\right)^{3} - 1}{3 \\left(z - \\frac{z^{3} - 1}{3 z^{2}}\\right)^{2}} - \\frac{z^{3} - 1}{3 z^{2}} - r_i\\,\\right|<\\epsilon,\\,\\forall r_i(\\,f(r_i)=0 )\\right\\}\$"
```
"""
basin(K::Define,j) = K.newt ? nrset(K.E,K.m,j) : jset(K.E,j)

# define function for computing orbit of a z0 input
function orbit(K::Define{FT,QT,CT,M,N},z0::Complex{Float64}) where {FT,QT,CT,M,N}
    M ? (z = K.seed) : (z = z0)
    zn = 0x00
    while (N ? (K.Q(z,z0)::Float64>K.ϵ)::Bool : (K.Q(z,z0)::Float64<K.ϵ))::Bool && K.N>zn
        z = K.F(z,z0)::Complex{Float64}
        zn+=0x01
    end; #end
    # return the normalized argument of z or iteration count
    return (zn::UInt8,z::Complex{Float64})
end

"""
    Compute(::Fatou.Define)::Union{Matrix{UInt8},Matrix{Float64}}

`Compute` the `Array` for `Fatou.FilledSet` as specefied by `Fatou.Define`.
"""
function Compute(K::Define{FT,QT,CT,M,N})::Tuple{Matrix{UInt8},Matrix{Complex{Float64}}} where {FT,QT,CT,M,N}
    # generate coordinate grid
    Kyn = round(UInt16,(K.∂[4]-K.∂[3])/(K.∂[2]-K.∂[1])*K.n)
    x = range(K.∂[1]+0.0001,stop=K.∂[2],length=K.n)
    y = range(K.∂[4],stop=K.∂[3],length=Kyn)
    Z = x' .+ im*y # apply Newton-Orbit function element-wise to coordinate grid
    (matU,matF) = (Array{UInt8,2}(undef,Kyn,K.n),Array{Complex{Float64},2}(undef,Kyn,K.n))
    @time @threads for j = 1:length(y); for k = 1:length(x);
        (matU[j,k],matF[j,k]) = orbit(K,Z[j,k])::Tuple{UInt8,Complex{Float64}}
    end; end
    return (matU,matF)
end

import PyPlot: plot

function plot(K::FilledSet;c::String="",bare::Bool=false)
    # plot figure using imshow based in input preferences
    figure()
    isempty(c) && (c = K.meta.cmap)
    isempty(c) ? imshow(K.meta.iter ? K.iter : K.mix, extent=K.meta.∂) :
        imshow(K.meta.iter ? K.iter : K.mix, cmap=c, extent=K.meta.∂)
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

println("Fatou detected $(Threads.nthreads()) julia threads.")

end # module
