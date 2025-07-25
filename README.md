# PrettyStreamlines.jl

An extension to the Plots.jl ecosystem providing evenly‑spaced streamline plots of arbitrary density.
Under the hood it implements the classic Jobard–Lefer algorithm. 

---

## Features

* **Grid‑based fields**: supply `u` and `v` as matrices on an `x×y` mesh
* **Functional fields**: supply `u(x,y)` and `v(x,y)` as functions; they will be sampled on your mesh
* Configurable **minimum** and **maximum** streamline density

---

## Installation

```julia
using Pkg
Pkg.add("EvenStreamlines")
```

---

## Usage

```julia
using Plots
using EvenStreamlines
```

### 1. Grid‑based fields

```julia
# define a mesh
x = range(-1,1,length=30)
y = range(-1,1,length=30)

# simple rotational field: U, V are ny×nx matrices
U = [ -yj for yj in y, xi in x ]
V = [  xi for yj in y, xi in x ]

# plot with defaults
streamlines(x, y, U, V)

# customize density and styling
streamlines(x, y, U, V;
    min_density = 4,
    max_density = 12,
    lw = 2,
    c = :blue)
```

### 3. Functional fields

```julia
# define continuous vector field
u_fun(x,y) = sin(pi*x)*cos(pi*y)
v_fun(x,y) = -cos(pi*x)*sin(pi*y)

# EvenStreamlines will sample these on your mesh
streamlines(x, y, u_fun, v_fun;
    min_density = 3,
    max_density = 8,
    lw = 1)
```

---


### References

\[1] Jobard, B., & Lefer, W. (1997). *Creating Evenly‑Spaced Streamlines of Arbitrary Density*. In W. Lefer & M. Grave (Eds.), Visualization in Scientific Computing ’97: Proceedings of the Eurographics Workshop in Boulogne‑sur‑Mer, France, April 28–30, 1997 (pp. 43–55). Vienna: Springer. [https://doi.org/10.1007/978-3-7091-6876-9\_5](https://doi.org/10.1007/978-3-7091-6876-9_5)

\[2] Ma, K. (2025). *Evenly Spaced Streamlines* ([https://github.com/keithfma/evenly\_spaced\_streamlines](https://github.com/keithfma/evenly_spaced_streamlines)). GitHub. Retrieved July 25, 2025.

\[3] keithfma. (2025). *evenly\_spaced\_streamlines* ([https://github.com/keithfma/evenly\_spaced\_streamlines](https://github.com/keithfma/evenly_spaced_streamlines)). GitHub. Retrieved July 25, 2025.
