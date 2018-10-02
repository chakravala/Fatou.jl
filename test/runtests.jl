using Reduce,Fatou
using SyntaxTree
using LaTeXStrings
using PyPlot
using Test

load_package(:rlfi)
@test orbit(newton(:(z^3-1))) == nothing
@test orbit(:(1*x),x->x, [-1.7 2π -1.2],17,3,147) == nothing
@test orbit(:(1*x),x->x, [-1.7 2π -1.2],0,1) == nothing
@test basin(juliafill(:(z^3-1)),1) |> typeof == LaTeXStrings.LaTeXString
@test basin(newton(:(z^3-1)),1) |> typeof == LaTeXStrings.LaTeXString
@test (f=newton(:(z^3-1))|>fatou; plot(f) |> typeof == PyPlot.PyObject)
@test (f=mandelbrot(:(z^2+c))|>fatou ; plot(f) |> typeof == PyPlot.PyObject)
@test (f=juliafill(:(z^2-0.06+0.67im))|>fatou ; plot(f) |> typeof == PyPlot.PyObject)
@test (η = Fatou.newton_raphson(:(z^3-1),1); ξ = genfun(η,[:z,:c]); ζ = 2.1; [ξ(ζ,0), [Fatou.recomp(η,ζ,k) for k ∈ 2:3]...]) |> typeof <: Array
# fix for case k ∈ 2:4+
