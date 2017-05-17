module Fatou
using SymPy,PyPlot

#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export julia, mandelbrot, newton, nrset, plot, fatou

funk(r::String) = :((z,c)->$(parse(r)))
funK(r::String) = :((z,n,e)->$(parse(r)))

abstract AbstractFatou
type FatouMeta <: AbstractFatou
  F::Function # primary map
  Q::Function # escape criterion
  C::Function # complex fixed point coloring
  ∂::Array{Number,1} # bounds
  n::Integer # number of grid points
  N::Integer # number of iterations
  ϵ::Float64 # epsilon Limit criterion
  m::Number # newton multiplicity factor
  e::Number # iteration color exponent
  iter::Bool # toggle iteration mode
  newt::Bool # toggle Newton mode
  mandel::Bool # toggle Mandelbrot mode
  orig::Function # original Newton map
  start # orbit starting point
  orbit::Int # orbit cobweb depth
  depth::Int # depth of function composition
  cmap::String
  function FatouMeta(F::String;
      Q::String="abs2(z)",
      C::String="angle(z)/(2π))*n^e",
      ∂=π/2,
      n::Integer=176,
      N::Integer=35,
      ϵ::Number=4,
      m::Number=0,
      e::Number=0,
      iter::Bool=false,
      newt::Bool=false,
      mandel::Bool=false,
      orig::String=F,
      start=nothing,
      orbit::Int=0,
      depth::Int=1,
      cmap::String="")
    typeof(∂) ≠ Array{Float64,1} && (∂ = [-float(∂),∂,-∂,∂])
    !newt ? (f = funk(F) |> eval; q = funk(Q) |> eval) :
      (f = newton_raphson(eval(funk(F)),m); q = eval(funk("abs("*F*")")))
    c = funK(C) |> eval; o = funk(orig) |> eval
    return new(f,q,c,∂,n,N,float(ϵ),m,e,iter,newt,mandel,o,start,orbit,depth,cmap);
  end; end

immutable FatouSet <: AbstractFatou
  meta::FatouMeta
  set::Union{Matrix{UInt8},Matrix{Float64}}
  FatouSet(K::FatouMeta) = new(K,FatouComp(K)); end

fatou(K::FatouMeta) = FatouSet(K)

function julia(F::String;
    Q::String= "abs2(z)",
    C::String= "(angle(z)/(2π))*n^e",
    ∂=π/2,
    n::Integer=176,
    N::Integer=35,
    ϵ::Number=4,
    m::Number=0,
    e::Number=0,
    iter::Bool=true,
    newt::Bool=false,
    start=nothing,
    orbit::Int=0,
    depth::Int=1,
    cmap::String="")
  return FatouMeta(F,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,m=m,e=e,iter=iter,newt=newt,start=start,orbit=orbit,depth=depth,cmap=cmap); end

function mandelbrot(F::String;
    Q::String= "abs2(z)",
    C::String= "exp(-abs(z))",
    ∂=π/2,
    n::Integer=176,
    N::Integer=35,
    ϵ::Number=4,
    m::Number=0,
    e::Number=0,
    iter::Bool=true,
    newt::Bool=false,
    start=nothing,
    orbit::Int=0,
    depth::Int=1,
    cmap::String="")
  return FatouMeta(F,Q=Q,C=C,∂=∂,n=n,N=N,ϵ=ϵ,m=m,e=e,iter=iter,newt=newt,mandel=true,start=start,orbit=orbit,depth=depth,cmap=cmap); end

function newton(F::String;
    C::String= "(angle(z)/(2π))*n^e",
    ∂=π/2,
    n::Integer=176,
    N::Integer=35,
    ϵ::Number=0.01,
    m::Number=1,
    e::Number=0,
    iter::Bool=true,
    mandel::Bool=false,
    start=nothing,
    orbit::Int=0,
    depth::Int=1,
    cmap::String="")
  return FatouMeta(F,C=C,∂=∂,n=n,N=N,ϵ=ϵ,m=m,e=e,iter=iter,newt=true,mandel=mandel,orig=F,start=start,orbit=orbit,depth=depth,cmap=cmap); end

# generate function code by constructing the lambda expression
sym2fun(expr,typ) = Expr(:function, Expr(:call, gensym(), map(s->Expr(:(::),s,typ),sort!(Symbol.(free_symbols(expr))))..., Expr(:(...),:zargs)), SymPy.walk_expression(expr))

# evaluate the expression and assign to a handle
p1(f::Function) = sym2fun(diff(f(Sym(:a),Sym(:b))),:Any) |> eval

# we can substitute the expression into Newton's method and display it with LaTeX
newton_raphson(f::Function,m) = (z,c) -> z - m*f(z,c)/p1(f)(z,c)

# define recursive composition on functions
recomp(f::Function,x::Number,j::Int) = j > 1 ? f(recomp(f,x,j-1),0) : f(x,0)

# we can convert the j-th function composition into a latex expresion
nL(f::Function,m,j) = recomp(newton_raphson(f,m),Sym(:z),j) |> SymPy.latex

# set of points that are within an ϵ neighborhood of the roots ri of the function f
setstr = "- r_i\\,\\right|<\\epsilon,\\,\\forall r_i(\\,f(r_i)=0 )\\right\\}"
latexstring("D_0(\\epsilon) = \\left\\{ z\\in\\mathbb{C}: \\left|\\,z $setstr")

# each subsequent iteration of the Newton method will yield a more complicated set
nrset(f::Function,m,j) = latexstring(
  "$ds D_$j(\\epsilon) = \\left\\{z\\in\\mathbb{C}:\\left|\\,$(nL(f,m,j)) $setstr")
nrset(K::FatouMeta,j) = nrset(K.F,K.m,j)

include("orbitplot.jl"); ds = "\\displaystyle"

function FatouComp(K::FatouMeta)
  # define Complex{Float64} versions of polynomial and constant for speed
  f = sym2fun(K.F(Sym(:a),Sym(:b)),:(Complex{Float64})) |> eval
  h = sym2fun(K.Q(Sym(:a),Sym(:b)),:(Complex{Float64})) |> eval
  # define function for computing orbit of a z0 input
  function nf(z0::Complex{Float64})::Number
    K.mandel ? (z = 0.0 + 0.0im): (z = z0); zn = 0
    while (K.newt ? (h(z,z0)>K.ϵ) : (h(z,z0)<K.ϵ)) && K.N>zn
      z = f(z,z0); zn+=1; end; #end
    # return the normalized argument of z or iteration count
    return Number(K.iter ? UInt8(zn) : K.C(z,zn/K.N,K.e)); end
  # generate coordinate grid
  x = linspace(K.∂[1]+0.0001,K.∂[2],K.n); y = linspace(K.∂[3],K.∂[4],K.n)
  # apply Newton-Orbit function element-wise to coordinate grid
  return @time nf.(x' .+ im*y); end

import PyPlot: plot

function plot(K::FatouSet;c::String="")
  # plot figure using imshow based in input preferences
  figure(); isempty(c) && (c = K.meta.cmap)
  isempty(c) ? imshow(K.set,extent=K.meta.∂) : imshow(K.set,cmap=c,extent=K.meta.∂)
  # determine if plot is Iteration, Roots, or Limit
  typeof(K.set) == Matrix{UInt8} ? t = L"iter. " :
    K.meta.m==1 ? t = L"roots" : t = L"limit"
    # annotate title using LaTeX
    if K.meta.newt
      title(latexstring("f:z\\mapsto $(SymPy.latex(K.meta.orig(Sym(:z),Sym(:c)))),\\, m = $(K.meta.m), ")*t)
      # annotate y-axis with Newton's method
      ylabel(L"Fatou\,set:\,"*L"z\,↦\,z-m\,×\,f(z)\,/\,f\,'(z)")
    else
      title(latexstring("f:z\\mapsto $(SymPy.latex(K.meta.F(Sym(:z),Sym(:c)))),\\,")*t)
    end
    colorbar(); tight_layout(); end

end # module
