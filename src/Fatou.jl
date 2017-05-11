module Fatou
using SymPy,PyPlot

#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

export newton, nrset, NewtonFractal, PlotNF

include("orbitplot.jl"); ds = "\\displaystyle"

# generate function code by constructing the lambda expression
sym2fun(expr,typ) = Expr(:function, Expr(:call, gensym(), Expr(:(::), map(Symbol,free_symbols(expr))...,typ)), SymPy.walk_expression(expr))

# evaluate the expression and assign to a handle
p1(f::Function) = eval(sym2fun(diff(f(Sym("z"))),:Any))

# we can substitute the expression into Newton's method and display it with LaTeX
newton(f::Function,m) = z -> z - m*f(z)/p1(f)(z)

# define multiplication operator on functions to evaluate as function composition
import Base: *; *(f::Function,g::Function) = x -> f(g(x))

# we can convert the j-th function composition into a latex expresion
nL(f::Function,m,j) = SymPy.latex(j==1 ?
  newton(f,m)(Sym("z")) : (newton(f,m)^j)(Sym("z")))

# set of points that are within an ϵ neighborhood of the roots ri of the function f
setstr = "- r_i\\,\\right|<\\epsilon,\\,\\forall r_i(\\,f(r_i)=0 )\\right\\}"
latexstring("D_0(\\epsilon) = \\left\\{ z\\in\\mathbb{C}: \\left|\\,z $setstr")

# each subsequent iteration of the Newton method will yield a more complicated set
nrset(f::Function,m,j) = latexstring(
  "$ds D_$j(\\epsilon) = \\left\\{z\\in\\mathbb{C}:\\left|\\,$(nL(f,m,j)) $setstr")

function NewtonFractal(m::Number, p::Function, γ::Number, n::Int=176; N::Int=35, ϵ::Float64=0.,exp::Number=0,iter::Bool=false,col::Function=(x, y, z)->x*y^z)
  return NewtonFractal(m,p,[-γ,γ,-γ,γ],n,N,ϵ,exp,iter,col); end

function NewtonFractal(m::Number, p::Function, ∂::Array{Float64,1}, n::Int=176; N::Int=35, ϵ::Float64=0.,exp::Number=0,iter::Bool=false,col::Function=(x, y, z)->x*y^z)
  # define Complex{Float64} versions of polynomial and constant for speed
  p0(z::Complex{Float64}) = p(z); m = convert(Complex{Float64},m)
  p1 = eval(sym2fun(diff(p(Sym("z"))),:(Complex{Float64})))
  newton(z::Complex{Float64})::Complex{Float64} = z - m*p0(z)/p1(z)
  # define function for computing orbit of a z0 input
  function nf(z::Complex{Float64})::Number
    if ϵ==0; for zn ∈ 1:N; z = newton(z); end # loop over orbit of z0
    else zn=0; while abs(p0(z))>ϵ && N>zn
      z = newton(z); zn+=1; end; end
    # return the normalized argument of z or iteration count
    return Number(iter ? zn : col(angle(z)/(2π),zn/N,exp)); end
  # generate coordinate grid
  x = linspace(∂[1]+0.0001,∂[2],n); y = linspace(∂[3],∂[4],n)
  # apply Newton-Orbit function element-wise to coordinate grid
  return @time nf.(x' .+ im*y); end

function PlotNF(nf::Union{Array{Float64,2},Array{Int,2}}, ∂=0, f::Function=0, m::Number=1; c::String="")
  # plot figure using imshow based in input preferences
  figure(); if ∂==0; isempty(c) ? imshow(nf) : imshow(nf,cmap=c)
    else typeof(∂)!=Array{Float64,1} && (∂=[-∂,∂,-∂,∂]);
      isempty(c) ? imshow(nf,extent=∂) : imshow(nf,cmap=c,extent=∂); end
    # determine if plot is Iteration, Roots, or Limit
    typeof(nf) == Matrix{Int64} ? t = L", iter. " :
      m==1 ? t = L", roots" : t = L", limit"
    # annotate title using LaTeX
    f!=0 && title(latexstring("z\\mapsto $(SymPy.latex(f(Sym("z")))),\\, m = $m")*t);
    # annotate y-axis with Newton's method
    ylabel(L"Fatou\,set:\,"*L"z\,↦\,z-m\,×\,f(z)\,/\,f\,'(z)")
    colorbar(); tight_layout(); end

end # module
