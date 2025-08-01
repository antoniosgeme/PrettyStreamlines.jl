using RecipesBase

"""
`streamlines(xs, ys, u, v; kwargs...)` — Plot evenly-spaced streamlines of a 2D vector field.

# Arguments
- `xs, ys`  
  1D coordinate vectors or 2D meshgrid matrices defining the domain.
- `u, v`  
  Velocity components: either Matrices of size `(Ny x Nx)` or functions `u(x,y)` and `v(x,y)`.

# Keywords
- `min_density::Real = 3`  
  Minimum streamline density (controls spacing of start points).
- `max_density::Real = 10`  
  Maximum streamline density.
- `glyphs::Bool = false`  
  If `true`, draw arrow glyphs along each streamline.
- `arrow_every::Int = 10`  
  Place an arrow on every Nth segment when `glyphs=true`.
- `arrow_scale::Real = 0.1`  
  Scale factor for arrow length.
- `color`  
  Fixed line color (e.g. `:black`), ignored if `color_by` is set.
- `color_by`  
  Either `:magnitude` or a function `f(x,y,u,v)` returning a value to map to color.
- `unbroken::Bool = false`  
  If `true`, do not truncate lines upon collision with existing streamlines.
- `seeds = nothing`  
  Optional seed points as a `Vector{Tuple{T,T}}` or an `N×2` matrix; if `nothing`, seeds are auto-generated.

# Examples

1. Basic usage with functions:
```jldoctest
julia> using Plots, PrettyStreamlines

julia> x = LinRange(-1,1,50)
julia> y = LinRange(-1,1,50)
julia> u(x,y) = -y  
julia> v(x,y) = x
julia> streamlines(x, y, u, v; min_density=2, max_density=8, color=:blue))
```                    
"""
@userplot Streamlines

@recipe function f(streams::Streamlines;
                  min_density   = 3,
                  max_density   = 10,
                  glyphs        = false,
                  arrow_every   = 10,
                  arrow_scale   = 0.1,
                  color         = :black,      # fixed color
                  color_by      = nothing,      # :magnitude or Function
                  unbroken      = false,
                  seeds         = nothing
)
    x, y, u, v = streams.args
    X,Y,U,V = process_stream_fields(x,y,u,v)

    # --- trace streamlines ---
    xy = compute_streamlines(x, y, u, v; min_density=min_density, 
                                     max_density=max_density,
                                     unbroken=unbroken,
                                     seeds=seeds)

    # --- prepare quiver arrows ---
    valid = .!isnan.(xy[:,1])
    idxs  = findall(valid)
    safe  = idxs[(arrow_every+1) : arrow_every : end-arrow_every]

    xs = Float64[]; ys = Float64[]
    dx = Float64[]; dy = Float64[]
    for j in safe
        x0,y0 = xy[j-1,1], xy[j-1,2]
        x1,y1 = xy[j+1,1], xy[j+1,2]
        tx,ty = x1-x0, y1-y0
        n = hypot(tx,ty)
        push!(xs, xy[j,1]); push!(ys, xy[j,2])
        if n == 0
            push!(dx, 0); push!(dy, 0)
        else
            push!(dx, tx/n*arrow_scale)
            push!(dy, ty/n*arrow_scale)
        end
    end

    # nearest‐neighbor lookup
    nearest_index(arr, v) = argmin(abs.(arr .- v))

    # --- decide if we’re color‐mapping ---
    has_map = false
    color_fn = nothing

    if color_by isa AbstractVector
        @assert length(color_by) == size(xy,1) "color_by array must match total points"
        has_map = true
        cvals   = color_by

    elseif color_by === :magnitude
        has_map = true
        cvals = similar(xy[:,1])
        for (i,(xi, yi)) in enumerate(eachrow(xy))
            if isnan(xi)
                cvals[i] = NaN
            else
                ix = argmin(abs.(x .- xi))
                iy = argmin(abs.(y .- yi))
                cvals[i] = hypot(U[iy,ix], V[iy,ix])
            end
        end

    elseif color_by isa Function
        has_map = true
        cvals = similar(xy[:,1])
        for (i,(xi, yi)) in enumerate(eachrow(xy))
            if isnan(xi)
                cvals[i] = NaN
            else
                ix = argmin(abs.(x .- xi))
                iy = argmin(abs.(y .- yi))
                cvals[i] = color_by(xi, yi, U[iy,ix], V[iy,ix])
            end
        end
    end


    if color_by === :magnitude
        has_map = true
        color_fn = (x,y,u,v)->hypot(u,v)

    elseif color_by isa Function
        has_map = true
        color_fn = color_by
    end

    # --- compute per‐point values if needed ---
    cvals = nothing
    if has_map
        c = similar(xy[:,1])
        for (i,(xi, yi)) in enumerate(eachrow(xy))
            if isnan(xi)
                c[i] = NaN
            else
                ix = nearest_index(x, xi)
                iy = nearest_index(y, yi)
                ui = U[iy, ix]
                vi = V[iy, ix]
                c[i] = color_fn(xi, yi, ui, vi)
            end
        end
        cvals = c
    end

    arrow_c = has_map ? cvals[safe] : nothing

    # --- main streamline series ---
    @series begin
        if has_map
            line_z := cvals
        else
            color  := color         # fixed
        end
        (xy[:,1], xy[:,2])
    end


    if glyphs
        @series begin
            seriestype := :quiver
            quiver     := (dx, dy)
            arrow      := :closed
            if has_map
                color := nothing
                line_z  := repeat(arrow_c, inner=4)
            else
                color := color
            end
            (xs, ys)
        end
    end
end
