module Fatou
using SyntaxTree,Reduce,LaTeXStrings,Requires,Base.Threads

#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export fatou, juliafill, mandelbrot, newton, basin, plot

abstract type ComplexBundle end

struct Rectangle <: ComplexBundle
    ∂::Vector{Float64} # bounds
    n::UInt16 # number of grid points
    function Rectangle(∂::T=π/2,n::Integer=176) where T
        !(T <: Array) && (∂ = [-float(∂),∂,-∂,∂])
        length(∂) == 2 && (∂ = [∂[1],∂[2],∂[1],∂[2]])
        new(convert(Vector{Float64},∂),UInt16(n))
    end
end

struct ComplexRectangle <: ComplexBundle
    ∂::Rectangle
    Ω::Matrix{Complex{Float64}}
end

bounds(R::ComplexRectangle) = bounds(Rectangle(R))
bounds(Ω::Rectangle) = Ω.∂
Base.size(Ω::Rectangle) = (round(UInt16,(Ω.∂[4]-Ω.∂[3])/(Ω.∂[2]-Ω.∂[1])*Ω.n),Ω.n)
Base.size(R::ComplexRectangle) = size(R.Ω)

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
      cmap::String= "",                   # imshow color map
      plane::Bool = false,                # convert input disk to half-plane
      disk::Bool  = false)                # convert output half-plane to disk

`Define` the metadata for a `Fatou.FilledSet`.
"""
struct Define{FT<:Function,QT<:Function,CT<:Function,M,N,P,D} <: ComplexBundle
    E::Any # input expression
    F::FT # primary map
    Q::QT # escape criterion
    C::CT # complex fixed point coloring
    Ω::Rectangle # bounded grid points
    N::UInt16 # number of iterations
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
    plane::Bool # convert input disk to half-plane
    disk::Bool # convert output half-plane to disk
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
            cmap::String="",
            plane::Bool=false,
            disk::Bool=false)
        !newt ? (f = genfun(E,[:z,:c]); q = genfun(Q,[:z,:c])) :
        (f = genfun(newton_raphson(E,m),[:z,:c]); q = genfun(Expr(:call,:abs,E),[:z,:c]))
        c = genfun(C,[:z,:n,:p])
        e = typeof(E) == String ? parse(E) : E
        return new{typeof(f),typeof(q),typeof(c),mandel,newt,plane,disk}(e,f,q,c,Rectangle(∂,n),UInt16(N),float(ϵ),iter,float(p),newt,m,mandel,seed,x0,orbit,depth,cmap,plane,disk)
    end
end

"""
    Fatou.FilledSet(::Fatou.Define)

Compute the `Fatou.FilledSet` set using `Fatou.Define`.
"""
struct FilledSet{FT,QT,CT,M,N,P,D} <: ComplexBundle
    meta::Define{FT,QT,CT,M,N,P,D}
    set::ComplexRectangle
    iter::Matrix{UInt16}
    mix::Matrix{Float64}
    function FilledSet{FT,QT,CT,M,N,P,D}(K::Define{FT,QT,CT,M,N,P,D},Z::ComplexRectangle) where {FT,QT,CT,M,N,P,D}
        (i,s) = Compute(K,Z)
        return new{FT,QT,CT,M,N,P,D}(K,s,i,broadcast(K.C,s.Ω,broadcast(float,i./K.N),K.p))
    end
end

Rectangle(K::Define) = K.Ω
Rectangle(K::FilledSet) = Rectangle(ComplexRectangle(K))
Rectangle(R::ComplexRectangle) = R.∂

ComplexRectangle(K::FilledSet) = K.set
ComplexRectangle(Ω::Matrix{Complex{Float64}}) = ComplexRectangle(Rectangle([0,size(Ω)[2],0,size(Ω)[1]],size(Ω)[2]),Ω)

bounds(K::Define) = bounds(Rectangle(K))
bounds(K::FilledSet) = bounds(ComplexRectangle(K))

ranges(K::FilledSet) = ranges(K.meta)
ranges(K::Define) = ranges(Rectangle(K))
function ranges(Ω::Rectangle)
    yn,xn = size(Ω)
    x = range(Ω.∂[1]+0.0001,stop=Ω.∂[2],length=xn)
    y = range(Ω.∂[4],stop=Ω.∂[3],length=yn)
    return x,y
end

(K::Define)(Z) = fatou(K,Z)
(K::FilledSet)(Z) = fatou(K,Z)

"""
      fatou(::Fatou.Define)

Compute the `Fatou.FilledSet` set using `Fatou.Define`.

# Examples
```Julia
julia> fatou(K)
```
"""
fatou(K::Define{FT,QT,CT,M,N,P,D},Z::ComplexRectangle) where {FT,QT,CT,M,N,P,D} = FilledSet{FT,QT,CT,M,N,P,D}(K,Z)
fatou(K::Define,Z::Rectangle=Rectangle(K)) = fatou(K,fatou(Z))
fatou(K::Define,Z::FilledSet) = fatou(K,ComplexRectangle(Z))
fatou(K::Define,Z::Define) = fatou(K,fatou(Z))
fatou(K::FilledSet,Z=ComplexRectangle(K)) = fatou(K.meta,Z)
function fatou(K::Rectangle) # generate coordinate grid
    x, y = ranges(K)
    ComplexRectangle(K, x' .+ im*y)
end # apply Newton-Orbit function element-wise to coordinate grid

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
      cmap::String= "",                   # imshow color map
      plane::Bool = false,                # convert input disk to half-plane
      disk::Bool  = false)                # convert output half-plane to disk

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
        cmap::String="",
        plane::Bool=false,
        disk::Bool=false)
    return Define(E,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=newt,m=m,x0=x0,orbit=orbit,depth=depth,cmap=cmap,plane=plane,disk=disk)
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
      cmap::String= "",                   # imshow color map
      plane::Bool = false,                # convert input disk to half-plane
      disk::Bool  = false)                # convert output half-plane to disk

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
        cmap::String="",
        plane::Bool=false,
        disk::Bool=false)
    m ≠ 0 && (newt = true)
    return Define(E,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=newt,m=m,mandel=true,seed=seed,x0=x0,orbit=orbit,depth=depth,cmap=cmap,plane=plane,disk=disk)
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
      cmap::String= "",                   # imshow color map
      plane::Bool = false,                # convert input disk to half-plane
      disk::Bool  = false)                # convert output half-plane to disk

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
        cmap::String="",
        plane::Bool=false,
        disk::Bool=false)
    return Define(E,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=true,m=m,mandel=mandel,seed=seed,x0=x0,orbit=orbit,depth=depth,cmap=cmap,plane=plane,disk=disk)
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

plane(z::Complex) = (2z.re/(z.re^2+(1-z.im)^2))+im*(1-z.re^2-z.im^2)/(z.re^2+(1-z.im)^2)
disk(z::Complex) = (2z.re/(z.re^2+(1+z.im)^2))+im*(z.re^2+z.im^2-1)/(z.re^2+(1+z.im)^2)

# define function for computing orbit of a z0 input
function orbit(K::Define{FT,QT,CT,M,N,P,D},z0::Complex{Float64}) where {FT,QT,CT,M,N,P,D}
    M ? (z = K.seed) : (z = P ? plane(z0) : z0)
    zn = 0x0000
    while (N ? (K.Q(z,z0)::Float64>K.ϵ)::Bool : (K.Q(z,z0)::Float64<K.ϵ))::Bool && K.N>zn
        z = K.F(z,z0)::Complex{Float64}
        zn+=0x0001
    end; #end
    # return the normalized argument of z or iteration count
    return (zn::UInt16,(D ? disk(z) : z)::Complex{Float64})
end

"""
    Compute(::Fatou.Define)::Union{Matrix{UInt8},Matrix{Float64}}

`Compute` the `Array` for `Fatou.FilledSet` as specefied by `Fatou.Define`.
"""
function Compute(K::Define{FT,QT,CT,M,N,D},Z::ComplexRectangle)::Tuple{Matrix{UInt16},ComplexRectangle} where {FT,QT,CT,M,N,D}
    yn,xn = size(Z)
    (matU,matF) = (Matrix{UInt16}(undef,yn,xn),Matrix{Complex{Float64}}(undef,yn,xn))
    @time @threads for j = 1:yn; for k = 1:xn;
        (matU[j,k],matF[j,k]) = orbit(K,Z.Ω[j,k])::Tuple{UInt16,Complex{Float64}}
    end; end
    return (matU,ComplexRectangle(Rectangle(Z),matF))
end

# determine if plot is Iteration, Roots, or Limit
typeplot(K::FilledSet) = typeof(K.meta.iter ? K.iter : K.mix) == Matrix{UInt16} ? "iter." : K.meta.m==1 ? "roots" : "limit"

function Base.String(K::FilledSet) # annotate title using LaTeX
    text,t = "f : z ↦ $(K.meta.E),",typeplot(K)
    K.meta.newt ? "$text m = $(K.meta.m), $t" : "$text $t"
end

function __init__()
    println("Fatou detected $(Threads.nthreads()) julia threads.")
    @require ColorSchemes="35d6a980-a343-548e-a6ea-1d62b119f2f4" begin
        nonan(x) = isnan(x) ? 0.0 : x
        function (C::ColorSchemes.ColorScheme)(K::FilledSet)
            S = size(K.iter)
            H = zeros(ColorSchemes.RGB{Float64},S...)
            if K.meta.iter
                M = length(C)/(maximum(K.iter)+1)
                for x ∈ 1:S[1], y ∈ 1:S[2]
                    H[x,y] = C[round(Int,M*(K.iter[x,y]+1),RoundUp)]
                end
            else
                for x ∈ 1:S[1], y ∈ 1:S[2]
                    H[x,y] = get(C,nonan(K.mix[x,y]))
                end
            end
            return H
        end
    end
    @require ImageInTerminal="d8c32880-2388-543b-8c61-d9f865259254" begin
        import ColorSchemes
        function Base.show(io::IO,K::FilledSet;c::String="",bare::Bool=false)
            isempty(c) && (c = K.meta.cmap)
            display(getproperty(ColorSchemes, isempty(c) ? :balance : Symbol(c))(K))
            !bare && print(io,String(K))
        end
    end
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("makie.jl")
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("pyplot.jl")
    @require UnicodePlots="b8865327-cd53-5732-bb35-84acbb429228" include("uniplots.jl")
end

end # module
