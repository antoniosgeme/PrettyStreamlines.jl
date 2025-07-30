
@inline function rect_exit(xmin::T, xmax::T, ymin::T, ymax::T,
                           x0::T, y0::T, x1::T, y1::T) where T<:Real
    dx = x1 - x0
    dy = y1 - y0

    # time to hit vertical sides
    t_x = if      dx > zero(T)
              (xmax - x0) / dx
          elseif  dx < zero(T)
              (xmin - x0) / dx
          else
              typemax(T)   # never hits a vertical
          end

    # time to hit horizontal sides
    t_y = if      dy > zero(T)
              (ymax - y0) / dy
          elseif  dy < zero(T)
              (ymin - y0) / dy
          else
              typemax(T)   # never hits a horizontal
          end

    # pick the smaller positive t
    t = t_x < t_y ? t_x : t_y

    return x0 + t*dx, y0 + t*dy
end


# Trace a single streamline using RK2
function trace(
    x0::Float64, y0::Float64,
    itp_u, itp_v,
    stepsize::Float64, maxvert::Int,
    xmin::Float64, xmax::Float64,
    ymin::Float64, ymax::Float64,
    sign::Int
)::Matrix{Float64}

    # Preallocate
    verts = Array{Float64}(undef, maxvert+1, 2)
    verts[1,1], verts[1,2] = x0, y0

    x, y = x0, y0
    n = 1

    @inbounds for i in 1:maxvert
        # --- midpoint RK2 step ---
        u1 = itp_u(y, x);  v1 = itp_v(y, x)

        if isnan(u1) || isnan(v1)
            break
        end

        xm = x + sign * u1 * stepsize/2
        ym = y + sign * v1 * stepsize/2

        if !(xmin <= xm <= xmax && ymin <= ym <= ymax)
            xi, yi = rect_exit(xmin, xmax, ymin, ymax, x, y, xm, ym)
            n += 1
            verts[n,1], verts[n,2] = xi, yi
            break

        elseif isnan(itp_u(ym, xm)) || isnan(itp_v(ym, xm))
            break
        end

        x_new = x + sign * itp_u(ym, xm) * stepsize
        y_new = y + sign * itp_v(ym, xm) * stepsize

        if !(xmin <= x_new <= xmax && ymin <= y_new <= ymax)
            xi, yi = rect_exit(xmin, xmax, ymin, ymax, x, y, x_new, y_new)
            n += 1
            verts[n,1], verts[n,2] = xi, yi
            break

        elseif isnan(itp_u(y_new, x_new)) || isnan(itp_v(y_new, x_new))
            break
        end

        # 5) accept the step
        x, y = x_new, y_new
        n += 1
        verts[n,1], verts[n,2] = x, y
    end

    return @view verts[1:n, :]
end

function streamRK2(
    X::AbstractMatrix{T}, Y::AbstractMatrix{T},
    U::AbstractMatrix{T}, V::AbstractMatrix{T},
    min_density::Real, max_density::Real,unbroken::Bool,
    seeds::Union{Vector{Tuple{T,T}},Nothing}) where T<:Real

    has_seeds = !isnothing(seeds)
    x_coords, y_coords = X[1, :], Y[:, 1]

    num     = 20
    nrstart = ceil(Int, num * min_density)
    ncstart = ceil(Int, num * min_density)
    nrend   = ceil(Int, num * max_density)
    ncend   = ceil(Int, num * max_density)

    xmin, xmax = minimum(X), maximum(X)
    ymin, ymax = minimum(Y), maximum(Y)
    xrange     = xmax - xmin
    yrange     = ymax - ymin

    incstartx = xrange / ncstart
    incstarty = yrange / nrstart
    ixrangecs = ncstart / xrange * (1 - eps())
    iyrangers = nrstart / yrange * (1 - eps())
    ixrangece = ncend   / xrange * (1 - eps())
    iyrangere = nrend   / yrange * (1 - eps())

    stepsize = min(xrange / (ncend * 2), yrange / (nrend * 2), 0.1) 
    maxvert  = min(10000, round(Int, sum(size(V)) * 4 / stepsize))

    startgrid = zeros(Bool, nrstart, ncstart)
    endgrid   = zeros(Bool, nrend,   ncend)
    rc_list   = collect(vec(CartesianIndices((1:nrstart,1:ncstart))))
    shuffle!(rc_list)

    seed_cells = CartesianIndex[]
    if has_seeds
        for seed in seeds
            xs, ys = seed
            r = clamp(floor(Int, (ys - ymin)/incstarty) + 1, 1, nrstart)
            c = clamp(floor(Int, (xs - xmin)/incstartx) + 1, 1, ncstart)
            push!(seed_cells, CartesianIndex(r,c))
        end
        seed_set = Set(seed_cells)
        seed_map = Dict(cell => i for (i,cell) in enumerate(seed_cells))
    else
        seed_set = rc_list
    end

    U[isinf.(U)] .= NaN 
    V[isinf.(V)] .= NaN

    itp_u = interpolate((y_coords, x_coords), U, Gridded(Linear()))
    itp_v = interpolate((y_coords, x_coords), V, Gridded(Linear()))

    function trim_path(v::Matrix{T}, unbroken::Bool=false) where T
        
        for i in 1:size(v,1)
            if isnan(v[i,1]) || isnan(v[i,2])
                v = @view v[1:i-1, :]
                break
            end
        end

        tcc = floor(Int, (v[1,1] - xmin) * ixrangece) + 1
        trr = floor(Int, (v[1,2] - ymin) * iyrangere) + 1

        for j in axes(v,1)
            xc, yc = v[j,1], v[j,2]

            cc_s = clamp(floor(Int, (xc - xmin)*ixrangecs) + 1, 1, ncstart)
            rr_s = clamp(floor(Int, (yc - ymin)*iyrangers) + 1, 1, nrstart)
            startgrid[rr_s, cc_s] = true

            cc_e = floor(Int, (xc - xmin)*ixrangece) + 1
            rr_e = floor(Int, (yc - ymin)*iyrangere) + 1

            if cc_e < 1 || cc_e > ncend || rr_e < 1 || rr_e > nrend
                return @view v[1:j, :]
            end

            if !unbroken && endgrid[rr_e, cc_e] && !(cc_e==tcc && rr_e==trr)
                return @view v[1:j, :]
            end

            endgrid[rr_e, cc_e] = true
            tcc, trr = cc_e, rr_e
        end
        return v
    end

    vertsout = Vector{Array{Float64,2}}()
    for idx in rc_list
        r, c = idx.I
        if !startgrid[r, c] && (idx in seed_set)
            startgrid[r, c] = 1
            if !has_seeds
                x0 = xmin + (c - 0.5) * incstartx
                y0 = ymin + (r - 0.5) * incstarty
            else
                i_seed = seed_map[idx]
                x0, y0 = seeds[i_seed]
            end

            vf = trace(x0, y0, itp_u, itp_v, stepsize, maxvert, xmin, xmax, ymin, ymax,  1)
            vb = trace(x0, y0, itp_u, itp_v, stepsize, maxvert, xmin, xmax, ymin, ymax, -1)

            if !unbroken
                vfo = trim_path(vf,false)
                vbo = trim_path(vb,false) 
                second_half = vbo[end:-1:1, :]
                first_half  = vfo[2:end, :]
                vo = vcat(second_half, first_half)
            else
                second_half = vb[end:-1:1, :]
                first_half  = vf[2:end, :]
                full_line = vcat(second_half, first_half)
                vo = trim_path(full_line,true)
            end 
            push!(vertsout, vo)
            push!(vertsout, [NaN NaN])

        end # if startgrid
    end # end idx
    if isempty(vertsout)
        return reshape(Float64[], 0, 2)
    else
        return reduce(vcat, vertsout)
    end
end 



"""
    get_streamlines(
        xx, yy, uu, vv;
        min_density::Real = 1.0,
        max_density::Real = 5.0,
        unbroken::Bool = false,
        seeds::Union{Nothing, Vector{Tuple{T,T}}, AbstractMatrix{T}} = nothing
    ) where T<:Real

Compute evenly‐spaced streamlines for the vector field `(uu, vv)` over the grid `(xx, yy)`.

# Arguments
- `xx, yy` : Either  
  - 1D coordinate vectors of length N (for x) and M (for y), or  
  - 2D matrices of size M×N defining a meshgrid.  
- `uu, vv` : Either  
  - M×N arrays of velocity components matching `(xx, yy)`, or  
  - Functions `uu(x,y)` and `vv(x,y)` that will be evaluated on the grid.

# Keyword Arguments
- `min_density::Real = 1.0`  
  Minimum streamline density (controls starting point spacing).
- `max_density::Real = 5.0`  
  Maximum streamline density.
- `unbroken::Bool = false`  
  If `true`, do not truncate streamlines upon collision with other lines.
- `seeds::Union{Nothing, Vector{Tuple{T,T}}, AbstractMatrix{T}} = nothing`  
  Optional user‐provided seed points.  
  - As a vector of `(x, y)` tuples, or  
  - As an N×2 matrix of coordinates.  
  If `nothing`, seeds are auto‑generated according to `min_density`.

# Returns
A 2‑column `Matrix{Float64}` of concatenated `(x, y)` points for all streamlines,
with `NaN` rows separating individual lines.

# Examples
```julia
# Using grid vectors and array‐defined field:
xs = LinRange(0, 1, 20)
ys = LinRange(0, 1, 20)
u = (x,y) -> -y
v = (x,y) -> x
xy = get_streamlines(xs, ys, u, v; min_density=1, max_density=5.0)

# Using explicit seed points:
seeds = [ (0.1, 0.1), (0.5, 0.5), (0.9, 0.2) ]
xy2 = get_streamlines(xs, ys, u, v; seeds=seeds, unbroken=true)
"""
function get_streamlines(xx, yy, uu, vv;
                        min_density::Real = 1.0,
                        max_density::Real = 5.0,
                        unbroken::Bool = false,
                        seeds::Union{Nothing, Vector{Tuple{T,T}}, AbstractMatrix{T}} = nothing
                        ) where T<:Real
    
    X,Y,U,V = process_stream_fields(xx,yy,uu,vv)
    @assert size(X) == size(Y) == size(U) == size(V)
    @assert min_density > 0 && max_density > 0

    seeds = normalize_seeds(seeds)

    xy = streamRK2(
        X, Y,      
        U, V,      
        min_density, max_density,
        unbroken,
        seeds
    )
    return xy
end 

normalize_seeds(::Nothing) = nothing

function normalize_seeds(seeds::Vector{Tuple{T,T}}) where T<:Real
    return seeds
end

function normalize_seeds(seeds::AbstractMatrix{T}) where T<:Real
    @assert size(seeds, 2) == 2 "if you pass a matrix for seeds, it must be N×2"
    return [(seeds[i,1], seeds[i,2]) for i in axes(seeds,1)]
end

function process_stream_fields(xx,yy,uu,vv)
    if isa(xx, AbstractVector) && isa(yy, AbstractVector)
        xc, yc = xx, yy
        X = [i for j in yc, i in xc]
        Y = [j for j in yc, i in xc]
    elseif isa(xx, AbstractMatrix) && isa(yy, AbstractMatrix)
        X, Y = xx, yy
    else
        throw(ArgumentError("Streamlines: xx,yy must be both vectors or both matrices"))
    end


    # 2) build U and V on that same grid
    if isa(uu, AbstractMatrix) && isa(vv, AbstractMatrix)
        U, V = uu, vv
    elseif  isa(uu, Function) && isa(vv, Function)
        U = [ uu(x, y) for y in Y[:, 1], x in X[1, :]]
        V = [ vv(x, y) for y in Y[:, 1], x in X[1, :]]
    else
        throw(ArgumentError("Streamlines: uu,vv must be matrices or functions"))
    end

    return X,Y,U,V
end 