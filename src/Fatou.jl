module Fatou
using SymPy,PyPlot,Base.Threads

#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export fatou, juliafill, mandelbrot, newton, basin, plot

abstract AbstractFatou

"""
    Fatou.Define(::String;                # primary map, (z, c) -> F
      Q::String   = "abs2(z)",            # escape criterion, (z, c) -> Q
      C::String   = "angle(z)/(2π))*n^p", # coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # vertical grid points
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
type Define <: AbstractFatou
  F::Function # primary map
  Q::Function # escape criterion
  C::Function # complex fixed point coloring
  ∂::Array{Float64,1} # bounds
  n::UInt16 # number of grid points
  N::UInt8 # number of iterations
  ϵ::Float64 # epsilon Limit criterion
  iter::Bool # toggle iteration mode
  p::Float64 # iteration color exponent
  newt::Bool # toggle Newton mode
  m::Number # newton multiplicity factor
  O::Function # original Newton map
  mandel::Bool # toggle Mandelbrot mode
  seed::Number # Mandelbrot seed value
  x0 # orbit starting point
  orbit::Int # orbit cobweb depth
  depth::Int # depth of function composition
  cmap::String # imshow color map
  function Define(F::String;
      Q::String="abs2(z)",
      C::String="angle(z)/(2π))*n^p",
      ∂=π/2,
      n::Integer=176,
      N::Integer=35,
      ϵ::Number=4,
      iter::Bool=false,
      p::Number=0,
      newt::Bool=false,
      m::Number=0,
      O::String=F,
      mandel::Bool=false,
      seed::Number=0.0+0.0im,
      x0=nothing,
      orbit::Int=0,
      depth::Int=1,
      cmap::String="")
    !(typeof(∂) <: Array) && (∂ = [-float(∂),∂,-∂,∂])
    length(∂) == 2 && (∂ = [∂[1],∂[2],∂[1],∂[2]])
    !newt ? (f = funk(F) |> eval; q = funk(Q) |> eval) :
      (f = newton_raphson(eval(funk(F)),m); q = eval(funk("abs("*F*")")))
    c = funK(C) |> eval; o = funk(O) |> eval
    return new(f,q,c,convert(Array{Float64,1},∂),UInt16(n),UInt8(N),float(ϵ),iter,float(p),newt,m,o,mandel,seed,x0,orbit,depth,cmap);
  end; end

"""
    Fatou.FilledSet(::Fatou.Define)

Compute the `Fatou.FilledSet` set using `Fatou.Define`.
"""
immutable FilledSet <: AbstractFatou
  meta::Define
  set::Union{Matrix{UInt8},Matrix{Float64}}
  FilledSet(K::Define) = new(K,Compute(K)); end

  """
      fatou(::Fatou.Define)

Compute the `Fatou.FilledSet` set using `Fatou.Define`.

# Examples
```Julia
julia> fatou(K)
```
"""
fatou(K::Define) = FilledSet(K)

"""
    juliafill(::String;                   # primary map, (z, c) -> F
      Q::String   = "abs2(z)",            # escape criterion, (z, c) -> Q
      C::String   = "angle(z)/(2π))*n^p", # coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # vertical grid points
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
function juliafill(F::String;
    Q::String= "abs2(z)",
    C::String= "(angle(z)/(2π))*n^p",
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
  return Define(F,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=newt,m=m,x0=x0,orbit=orbit,depth=depth,cmap=cmap); end

"""
    mandelbrot(::String;                  # primary map, (z, c) -> F
      Q::String   = "abs2(z)",            # escape criterion, (z, c) -> Q
      C::String   = "exp(-abs(z))*n^p",   # coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # vertical grid points
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
mandelbrot("z^2+c",n=800,N=20,∂=[-1.91,0.51,-1.21,1.21],cmap="nipy_spectral")
```
"""
function mandelbrot(F::String;
    Q::String= "abs2(z)",
    ∂=π/2,
    C::String= "exp(-abs(z))*n^p",
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
  return Define(F,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=newt,m=m,mandel=true,seed=seed,x0=x0,orbit=orbit,depth=depth,cmap=cmap); end

"""
    newton(::String;                      # primary map, (z, c) -> F
      C::String   = "angle(z)/(2π))*n^p", # coloring, (z, n=iter., p=exp.) -> C
      ∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
      n::Integer  = 176,                  # vertical grid points
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
function newton(F::String;
    C::String= "(angle(z)/(2π))*n^p",
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
  return Define(F,C=C,∂=∂,n=n,N=N,ϵ=ϵ,iter=iter,p=p,newt=true,m=m,O=F,mandel=mandel,seed=seed,x0=x0,orbit=orbit,depth=depth,cmap=cmap); end

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
basin(K::Define,j) = K.newt ? nrset(K.O,K.m,j) : jset(K.F,j)

"""
    Compute(::Fatou.Define)::Union{Matrix{UInt8},Matrix{Float64}}

`Compute` the `Array` for `Fatou.FilledSet` as specefied by `Fatou.Define`.
"""
function Compute(K::Define)::Union{Matrix{UInt8},Matrix{Float64}}
  # define Complex{Float64} versions of polynomial and constant for speed
  f = (sym2fun(K.F(Sym(:a),Sym(:b)),:(Complex{Float64})) |> eval)::Function
  h = (sym2fun(K.Q(Sym(:a),Sym(:b)),:(Complex{Float64})) |> eval)::Function
  c(z::Complex{Float64},n::Number,p::Number) = K.C(z,n,p)::Float64
  # define function for computing orbit of a z0 input
  function nf(z0::Complex{Float64})::Union{UInt8,Float64}
    K.mandel ? (z = K.seed): (z = z0); zn = 0x00
    while (K.newt ? (h(z,z0)::Float64>K.ϵ)::Bool : (h(z,z0)::Float64<K.ϵ))::Bool && K.N>zn
      z = f(z,z0)::Complex{Float64}; zn+=0x01; end; #end
    # return the normalized argument of z or iteration count
    return (K.iter ? zn::UInt8:c(z,float(zn/K.N),K.p)::Float64)::Union{UInt8,Float64}; end
  # generate coordinate grid
  Kyn = round(UInt16,(K.∂[4]-K.∂[3])/(K.∂[2]-K.∂[1])*K.n)
  x = linspace(K.∂[1]+0.0001,K.∂[2],K.n); y = linspace(K.∂[4],K.∂[3],Kyn)
  # apply Newton-Orbit function element-wise to coordinate grid
  mat = (K.iter ? Array{UInt8,2}(Kyn,K.n) : Array{Float64,2}(Kyn,K.n))
  @time @threads for j = 1:length(y); for k = 1:length(x);
      mat[j,k] = nf(x[k] + im*y[j]); end; end; return mat; end # nf.(x' .+ im*y)

import PyPlot: plot

function plot(K::FilledSet;c::String="",bare::Bool=false)
  # plot figure using imshow based in input preferences
  figure(); isempty(c) && (c = K.meta.cmap)
  isempty(c) ? imshow(K.set,extent=K.meta.∂) : imshow(K.set,cmap=c,extent=K.meta.∂)
  tight_layout(); if !bare
    # determine if plot is Iteration, Roots, or Limit
    typeof(K.set) == Matrix{UInt8} ? t = L"iter. " :
      K.meta.m==1 ? t = L"roots" : t = L"limit"
    # annotate title using LaTeX
    ttext = "f:z\\mapsto $(SymPy.latex(K.meta.O(Sym(:z),Sym(:c)))),\\,"
    if K.meta.newt; title(latexstring("$ttext m = $(K.meta.m), ")*t)
      # annotate y-axis with Newton's method
      ylabel(L"Fatou\,set:\,"*L"z\,↦\,z-m\,×\,f(z)\,/\,f\,'(z)")
    else; title(latexstring("$ttext")*t)
    end; tight_layout(); colorbar(); end; end

end # module
