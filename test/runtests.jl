using Fatou
using Base.Test

@test orbitplot(x->x, [-1.7 2π -1.2],17,3,147) == nothing
@test orbitplot(x->x, [-1.7 2π -1.2],0,1) == nothing
@test (η = newton(z->z,1); ζ = 2.1; [η(ζ), [(η^k)(ζ) for k ∈ 2:4]...]) |> typeof <: Array
@test nrset(z->z^3-1,1,1) == nrset(x->x^3-1,1,1)
@test (f=z -> z^3-1; m = 1; γ = π/2; nf = NewtonFractal(m,f,[-γ,γ,-γ,γ],200,ϵ=0.01)) |> typeof <: Array && PlotNF(nf,γ,f,m) == PlotNF(nf,γ,f,m)
