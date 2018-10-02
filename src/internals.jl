#   This file is part of Fatou.jl. It is licensed under the MIT license
#   Copyright (C) 2017 Michael Reed

import Base: invokelatest

rdpm(tex) = split(split(tex,"\n\\end{displaymath}")[1],"\\begin{displaymath}\n")[2]

# we can substitute the expression into Newton's method and display it with LaTeX
function newton_raphson(F,m)
    f = RExpr(F)
    return Algebra.:-(R"z",Algebra.:*(m,Algebra.:/(f,Algebra.df(f,:z)))) |> factor |> parse
end

# define recursive composition on functions
recomp(E,x,j::Int) = Algebra.sub((Expr(:(=),:z,j > 1 ? recomp(E,x,j-1) : x),:(c=0)),E)

# we can convert the j-th function composition into a latex expresion
nL(E,m,j) = recomp(newton_raphson(E,m),:z,j) |> Algebra.latex |> rdpm
jL(E,j) = recomp(E,:z,j) |> Algebra.latex |> rdpm

# set of points that are within an ϵ neighborhood of the roots ri of the function f
ds = "\\displaystyle"
set0 = "D_0(\\epsilon) = \\left\\{ z\\in\\mathbb{C}: \\left|\\,z"
setj(j) = "$ds D_$j(\\epsilon) = \\left\\{z\\in\\mathbb{C}:\\left|\\,"
nsetstr = "- r_i\\,\\right|<\\epsilon,\\,\\forall r_i(\\,f(r_i)=0 )\\right\\}"
jsetstr = "\\,\\right|>\\epsilon\\right\\}"
nset0 = latexstring("$set0 $nsetstr")
jset0 = latexstring("$set0 $jsetstr")
nrset(f,m,j) = latexstring(
    j == 0 ? "$set0 $nsetstr" : "$(setj(j))$(nL(f,m,j)) $nsetstr")
jset(f,j) = latexstring(
    j == 0 ? "$set0 $jsetstr" : "$(setj(j))$(jL(f,j)) $jsetstr")
