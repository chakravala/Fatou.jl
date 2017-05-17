using Fatou
using Base.Test

@test orbit(newton("z^3-1")) == nothing
@test orbit(x->x, [-1.7 2π -1.2],17,3,147) == nothing
@test orbit(x->x, [-1.7 2π -1.2],0,1) == nothing
@test (η = Fatou.newton_raphson((z,c)->z,1); ζ = 2.1; [η(ζ,0), [Fatou.recomp(η,ζ,k) for k ∈ 2:4]...]) |> typeof <: Array
@test nrset((z,c)->z^3-1,1,1) == nrset((x,c)->x^3-1,1,1)
@test (f=newton("z^3-1")|>fatou ; f|>typeof == Fatou.FatouSet && plot(f) == plot(f))
@test (f=mandelbrot("z^2+c")|>fatou ; f|>typeof == Fatou.FatouSet && plot(f) == plot(f))
@test (f=julia("z^2-0.06+0.67im")|>fatou ; f|>typeof == Fatou.FatouSet && plot(f) == plot(f))
