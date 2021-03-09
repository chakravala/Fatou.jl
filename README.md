# Fatou.jl

[![Build Status](https://travis-ci.org/chakravala/Fatou.jl.svg?branch=master)](https://travis-ci.org/chakravala/Fatou.jl) [![Build status](https://ci.appveyor.com/api/projects/status/mdathjmu7jg57u77?svg=true)](https://ci.appveyor.com/project/chakravala/fatou-jl) [![Coverage Status](https://coveralls.io/repos/github/chakravala/Fatou.jl/badge.svg?branch=master)](https://coveralls.io/github/chakravala/Fatou.jl?branch=master) [![codecov.io](http://codecov.io/github/chakravala/Fatou.jl/coverage.svg?branch=master)](http://codecov.io/github/chakravala/Fatou.jl?branch=master)

Julia package for Fatou sets. Install using `Pkg.add("Fatou")` in Julia. See [Explore Fatou sets & Fractals](https://crucialflow.com/Fatou.jl) in Wiki for detailed *examples*. This package provides: `fatou`, `juliafill`, `mandelbrot`, `newton`, `basin`, `plot`, and `orbit`; along with various internal functionality using `Reduce` and Julia expressions to help compute `Fatou.FilledSet` efficiently. Full documentation is included. The `fatou` function can be applied to a `Fatou.Define` object to produce a `Fatou.FilledSet`, which can then be passed as an argument to `plot` functions of `Makie`, `PyPlot`, `ImageInTerminal`. Creation of `Fatou.Define` objects is done via passing a `parse`-able function expression string (in variables `z`, `c`) and optional keyword arguments to `juliafill`, `mandelbrot`, and `newton`.

## Background

This package enables users of Julia lang to easily generate, explore, and share fractals of Julia, Mandelbrot, and Newton type. The name Fatou comes from the mathematician after whom the Fatou sets are named. Note that the Julia language is not named after the mathematician Julia after whom the Julia sets are named. This is a mere coincidence.

**Definition** *(Julia set)*: For any holomorphic function on a complex plane, the boundary of the set of points whose result diverges when the function is iteratively evaluated at each point.

**Definition** *(Fatou set)*: The Julia set’s complement is the set of fixed limit points from holomorphic recursion.

**Definition** *(Mandelbrot set)*: The set of points on a complex parameter space for which the holomorphic recursion does not go to infinity from a common starting point `z0`.

**Definition** *(Newton fractal)*: The Julia/Fatou set obtained from the recursion of the Newton method `z↦z−m⋅f(z)/f′(z)` applied to a holomorphic function.

The package has essentially two different plotting modes controlled by the `iter` boolean keyword, which toggles whether to color the image by iteration count or whether to use a default (or custom) limit-value coloring function.

The number of Julia threads available is detected at the startup and is reported it back. When a specified Fatou set is computed, multi-threading is used to compute the pixels.
Since each pixel is independent of any other pixel, it doesn’t matter in what order or on how many threads it is computed, the more you use the faster it is.
The environment variable `JULIA_NUM_THREADS` can be used to enable the multi-threading for more than 1 thread.

Please share with us your favorite fractals as `Fatou` code snippets!

## Examples

`Fatou.Define` provides the following optional keyword arguments:

```Julia
Q::Expr 	= :(abs2(z)),           # escape criterion, (z, c) -> Q
C::Expr 	= :((angle(z)/(2π))*n^p)# coloring, (z, n=iter., p=exp.) -> C
∂    = π/2, # Array{Float64,1}      # Bounds, [x(a),x(b),y(a),y(b)]
n::Integer  = 176,                  # vertical grid points
N::Integer  = 35,                   # max. iterations
ϵ::Number   = 4,                    # basin ϵ-Limit criterion
iter::Bool  = false,                # toggle iteration mode
p::Number   = 0,                    # iteration color exponent
newt::Bool  = false,                # toggle Newton mode
m::Number   = 0,                    # Newton multiplicity factor
mandel::Bool= false,                # toggle Mandelbrot mode
seed::Number= 0.0+0.0im,            # Mandelbrot seed value
x0          = nothing,              # orbit starting point
orbit::Int  = 0,                    # orbit cobweb depth
depth::Int  = 1,                    # depth of function composition
cmap::String= ""                    # imshow color map
```

A Fatou set is a collection of complex valued orbits of an iterated function. To help illustrate this, an additional feature is a plot function designed to visualize real-valued-orbits.
The program can be initialized with `using Fatou, PyPlot` or `Makie` or `ImageInTerminal`.
For `PyPlot` the `imshow` and `plot` methods can be used, while for `Makie` the `heatmap`, `contour`, `surface`, and `arrows` methods can be used.

When `using ImageInTerminal`, the display of a `Fatou.FilledSet` will be plotted automatically in the terminal.
The `orbit` method also has optional `UnicodePlots` compatibility.
Additional plotting support can be added via Pull-Request by adding another `Requires` script to the `__init__()` function definition.

The following is a cobweb orbit plot of a function:

```Julia
juliafill(:(z^2-0.67),∂=[-1.25,1.5],x0=1.25,orbit=17,depth=3,n=147) |> orbit
```

![img/orbit.png](img/orbit.png)

With `fatou` and `plot` it is simple to display a filled in Julia set:

```Julia
c = -0.06 + 0.67im
nf = juliafill(:(z^2+$c),∂=[-1.5,1.5,-1,1],N=80,n=1501,cmap="gnuplot",iter=true)
plot(fatou(nf), bare=true)
```

![img/filled-julia.png](img/filled-julia.png)

It is also possible to switch to `mandelbrot` mode:

```Julia
mandelbrot(:(z^2+c),n=800,N=20,∂=[-1.91,0.51,-1.21,1.21],cmap="gist_earth") |> fatou |> plot
```

![img/mandelbrot.png](img/mandelbrot.png)

`Fatou` also provides `basin` to display the the Newton / Fatou basins using set notation in LaTeX in `IJulia`.

```Julia
map(display,[basin(newton(:(z^3-1)),i) for i ∈ 1:3])
```

![D1(ϵ)](http://latex.codecogs.com/svg.latex?D_1(\epsilon)%20=%20\left\\{z\in\mathbb{C}:\left|\\,z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}%20-%20r_i\\,\right|%3C\epsilon,\\,\forall%20r_i(\\,f(r_i)=0%20)\right\\})

![D2(ϵ)](http://latex.codecogs.com/svg.latex?D_2(\epsilon)%20=%20\left\\{z\in\mathbb{C}:\left|\\,z%20-%20\frac{\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{3}%20-%201}{3%20\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{2}}%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}%20-%20r_i\\,\right|%3C\epsilon,\\,\forall%20r_i(\\,f(r_i)=0%20)\right\\})

![D3(ϵ)](http://latex.codecogs.com/svg.latex?D_3(\epsilon)%20=%20\left\\{z\in\mathbb{C}:\left|\\,z%20-%20\frac{\left(z%20-%20\frac{\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{3}%20-%201}{3%20\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{2}}%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{3}%20-%201}{3%20\left(z%20-%20\frac{\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{3}%20-%201}{3%20\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{2}}%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{2}}%20-%20\frac{\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{3}%20-%201}{3%20\left(z%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}\right)^{2}}%20-%20\frac{z^{3}%20-%201}{3%20z^{2}}%20-%20r_i\\,\right|%3C\epsilon,\\,\forall%20r_i(\\,f(r_i)=0%20)\right\\})

Compute the Newton fractal Julia set for a function with annotated plot of iteration count:

```Julia
nf = newton(:(z^3-1),n=800,ϵ=0.1,N=25,iter=true,cmap="jet")
nf |> fatou |> plot
basin(nf,3)
```

![img/newton.png](img/newton.png)

Generalized Newton fractal example:

```Julia
nf = newton(:(sin(z)-1),m=1-1im,∂=[-2π/3,-π/3,-π/6,π/6],n=500,N=33,iter=true,ϵ=0.05,cmap="cubehelix")
nf |> fatou |> plot
basin(nf,2)
```

![img/generalized-newton.png](img/generalized-newton.png)

![D2(ϵ)](http://latex.codecogs.com/svg.latex?D_2(\epsilon)%20=%20\left\\{z\in\mathbb{C}:\left|\\,z%20-%20\frac{1}{\cos{\left%20(z%20\right%20)}}%20\left(1%20+%20i\right)%20\left(\sin{\left%20(z%20\right%20)}%20-%201\right)%20-%20\frac{\left(1%20+%20i\right)%20\left(\sin{\left%20(z%20-%20\frac{1}{\cos{\left%20(z%20\right%20)}}%20\left(1%20+%20i\right)%20\left(\sin{\left%20(z%20\right%20)}%20-%201\right)%20\right%20)}%20-%201\right)}{\cos{\left%20(z%20-%20\frac{1}{\cos{\left%20(z%20\right%20)}}%20\left(1%20+%20i\right)%20\left(\sin{\left%20(z%20\right%20)}%20-%201\right)%20\right%20)}}%20-%20r_i\\,\right|%3C\epsilon,\\,\forall%20r_i(\\,f(r_i)=0%20)\right\\})

View [Explore Fatou sets & Fractals](https://crucialflow.com/Fatou.jl) in Wiki for detailed *examples*.

### Troubleshooting on Julia 1.0.1+

Note that `Fatou` is not compatible with Julia 1.0 but works on Julia 1.0.1 alright. Note that a stackoverflow error occurs on Julia 1.0.1+ when the `Reduce` package is precompiled with `ENV["REDPRE"]` flag set, therefore it is recommended to not set it.
If you encounter an unsatisfiable requirement in the package manager, an easy workaround is to use `dev Fatou` instead of `add Fatou`.
