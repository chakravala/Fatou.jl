module Fatou
using SymPy,PyPlot

#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export newton, nrset, NewtonFractal, PlotNF

include("orbitplot.jl"); ds = "\\displaystyle"

abstract AbstractFatou
type FatouMeta <: AbstractFatou
  F::Function # map
  L::Function
  C::Function
  ∂::Array{Number,1} # bounds
  n::Integer # number of points
  N::Integer # number of iterations
  ϵ::Float64 # epsilon Limit
  m::Number # newton multiplicity
  e::Number
  iter::Bool
  newt::Bool
  mandel::Bool
  function FatouMeta(F::Function;
      L::Function=(z,c)->abs2(z),
      ∂=π/2,
      n::Integer=176,
      N::Integer=35,
      ϵ::Number=4,
      m::Number=1,
      e::Number=0,
      C::Function=(z, n, e)->(angle(z)/(2π))*n^e,
      iter::Bool=false,
      newt::Bool=false,
      mandel::Bool=false)
    typeof(∂) ≠ Array{Float64,1} && (∂ = [-float(∂),∂,-∂,∂])
    return new(F,L,C,∂,n,N,float(ϵ),m,e,iter,newt,mandel)
  end
end
immutable FatouSet <: AbstractFatou
  meta::FatouMeta
  set::Union{Matrix{Integer},Matrix{Float64}}
  FatouSet(K::FatouMeta) = new(K,FatouComp(K))
end

function mandelbrot(F::Function,∂::Number;n::Integer=176,N::Integer=35)
  return FatouMeta(F,C=(z,n,e)->exp(-abs(z)),∂=∂,n=n,N=N,mandel=true)
end
#newton(), newton_raphson() julia()

function FatouComp(K::FatouMeta)
  # define Complex{Float64} versions of polynomial and constant for speed
  e = complex(float(K.e))
  f = sym2fun(K.F(Sym(:z),Sym(:c)),:(Complex{Float64})) |> eval
  h = sym2fun(K.L(Sym(:z),Sym(:c)),:(Complex{Float64})) |> eval
  c = sym2fun(K.C(Sym(:z),Sym(:n),Sym(:e)),:(Complex{Float64})) |> eval
  # define function for computing orbit of a z0 input
  function nf(z0::Complex{Float64})::Number
    K.mandel ? (z = 0.0im) : (z = z0)
    if K.ϵ==0; for zn ∈ 1:K.N; z = f(z,z0); end # loop over orbit of z0
    else zn=0; while h(z,z0)<K.ϵ && K.N>zn
      z = f(z,z0); zn+=1; end; end
    # return the normalized argument of z or iteration count
    return Number(K.iter ? zn : c(z,complex(float(zn/K.N)),e)); end
  # generate coordinate grid
  #if K.mandel; x=zeros(K.n); y = zeros(K.n)
  x = linspace(K.∂[1]+0.0001,K.∂[2],K.n); y = linspace(K.∂[3],K.∂[4],K.n)
  # apply Newton-Orbit function element-wise to coordinate grid
  return @time nf.(x' .+ im*y); end


# generate function code by constructing the lambda expression
sym2fun(expr,typ) = Expr(:function, Expr(:call, gensym(), map(s->Expr(:(::),Symbol(s),typ),free_symbols(expr))..., Expr(:(...),:args)), SymPy.walk_expression(expr))

# evaluate the expression and assign to a handle
p1(f::Function) = sym2fun(diff(f(Sym(:z))),:Any) |> eval

# we can substitute the expression into Newton's method and display it with LaTeX
newton_raphson(f::Function,m) = z -> z - m*f(z)/p1(f)(z)

# define recursive composition on functions
recomp(f::Function,x::Number,j::Int) = j > 1 ? f(recomp(f,x,j-1)) : f(x)

# we can convert the j-th function composition into a latex expresion
nL(f::Function,m,j) = recomp(newton_raphson(f,m),Sym(:z),j) |> SymPy.latex

# set of points that are within an ϵ neighborhood of the roots ri of the function f
setstr = "- r_i\\,\\right|<\\epsilon,\\,\\forall r_i(\\,f(r_i)=0 )\\right\\}"
latexstring("D_0(\\epsilon) = \\left\\{ z\\in\\mathbb{C}: \\left|\\,z $setstr")

# each subsequent iteration of the Newton method will yield a more complicated set
nrset(f::Function,m,j) = latexstring(
  "$ds D_$j(\\epsilon) = \\left\\{z\\in\\mathbb{C}:\\left|\\,$(nL(f,m,j)) $setstr")

function NewtonFractal(m::Number, p::Function, ∂, n::Int=176; N::Int=35, ϵ::Float64=0.,exp::Number=0,iter::Bool=false,col::Function=(x, y, z)->x*y^z)
  typeof(∂) ≠ Array{Float64,1} && (∂ = [-∂,∂,-∂,∂])
  # define Complex{Float64} versions of polynomial and constant for speed
  p0(z::Complex{Float64}) = p(z); m = convert(Complex{Float64},m)
  p1 = sym2fun(diff(p(Sym(:z))),:(Complex{Float64})) |> eval
  newton_raphson(z::Complex{Float64})::Complex{Float64} = z - m*p0(z)/p1(z)
  # define function for computing orbit of a z0 input
  function nf(z::Complex{Float64})::Number
    if ϵ==0; for zn ∈ 1:N; z = newton_raphson(z); end # loop over orbit of z0
    else zn=0; while abs(p0(z))>ϵ && N>zn
      z = newton_raphson(z); zn+=1; end; end
    # return the normalized argument of z or iteration count
    return Number(iter ? zn : col(angle(z)/(2π),zn/N,exp)); end
  # generate coordinate grid
  x = linspace(∂[1]+0.0001,∂[2],n); y = linspace(∂[3],∂[4],n)
  # apply Newton-Orbit function element-wise to coordinate grid
  return @time nf.(x' .+ im*y); end

function PlotNF(nf::Union{Array{Float64,2},Array{Int,2}}, ∂=0, f::Function=0, m::Number=1; c::String="")
  # plot figure using imshow based in input preferences
  figure(); if ∂==0; isempty(c) ? imshow(nf) : imshow(nf,cmap=c)
    else typeof(∂) ≠ Array{Float64,1} && (∂=[-∂,∂,-∂,∂])
      isempty(c) ? imshow(nf,extent=∂) : imshow(nf,cmap=c,extent=∂); end
    # determine if plot is Iteration, Roots, or Limit
    typeof(nf) == Matrix{Int64} ? t = L", iter. " :
      m==1 ? t = L", roots" : t = L", limit"
    # annotate title using LaTeX
    f≠0 && title(latexstring("f:z\\mapsto $(SymPy.latex(f(Sym(:z)))),\\, m = $m")*t);
    # annotate y-axis with Newton's method
    ylabel(L"Fatou\,set:\,"*L"z\,↦\,z-m\,×\,f(z)\,/\,f\,'(z)")
    colorbar(); tight_layout(); end

end # module
