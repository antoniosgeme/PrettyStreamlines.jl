using Interpolations: interpolate, Gridded, Linear
using Random: shuffle!

function meshgrid(x::AbstractVector, y::AbstractVector)
    X = [i for j in y, i in x]
    Y = [j for j in y, i in x]
    return X, Y
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
    # allocate the worst‑case number of points
    verts = Array{Float64}(undef, maxvert+1, 2)
    verts[1,1], verts[1,2] = x0, y0
    x, y = x0, y0
    n = 1
    @inbounds for i in 1:maxvert
        # exit early if out of bounds
        if !(xmin <= x <= xmax && ymin <= y <= ymax)
            break
        end
        # midpoint RK2
        u1 = itp_u(y, x);  v1 = itp_v(y, x)
        xm = x + sign * u1 * stepsize/2
        ym = y + sign * v1 * stepsize/2
        if !(xmin <= xm <= xmax && ymin <= ym <= ymax)
            break
        end
        u2 = itp_u(ym, xm);  v2 = itp_v(ym, xm)
        x += sign * u2 * stepsize
        y += sign * v2 * stepsize
        n += 1
        verts[n,1], verts[n,2] = x, y
    end
    return @view verts[1:n, :]  # no copy, returns a Matrix{Float64}
end

function get_stream_xy(
    X::AbstractMatrix{T}, Y::AbstractMatrix{T},
    U::AbstractMatrix{T}, V::AbstractMatrix{T},
    x_coords::AbstractVector{T}, y_coords::AbstractVector{T},
    min_density::Real, max_density::Real ) where T<:Real

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

    stepsize = min(xrange / (ncend * 2), yrange / (nrend * 2), 0.1) # A more robust stepsize
    maxvert  = min(10000, round(Int, sum(size(V)) * 4 / stepsize))

    startgrid = zeros(Bool, nrstart, ncstart)
    endgrid   = zeros(Bool, nrend,   ncend)
    rc_list   = collect(vec(CartesianIndices((1:nrstart,1:ncstart))))
    shuffle!(rc_list)

    U[isinf.(U)] .= NaN 
    V[isinf.(V)] .= NaN

    # Continuous interpolators
    itp_u = interpolate((y_coords, x_coords), U, Gridded(Linear()))
    itp_v = interpolate((y_coords, x_coords), V, Gridded(Linear()))

    function process_trace(v::Matrix{T})
        nan_idx = findfirst(isnan, v)
        if nan_idx !== nothing
            first_nan_row = nan_idx[1]
            if first_nan_row == 1
                v = v[1:0, :]
            else
                v = v[1:first_nan_row - 1, :]
            end
        end
        tcc = floor(Int, (v[1,1]-xmin)*ixrangece )+1
        trr = floor(Int, (v[1,2]-ymin)*iyrangere )+1
        break_j = size(v,1)
        for j=1:size(v,1)
            xc = v[j,1]
            yc = v[j,2]
            cc = floor(Int, (xc-xmin)*ixrangecs )+1
            rr = floor(Int, (yc-ymin)*iyrangers )+1
            if cc > 0 && cc <= ncstart && rr > 0 && rr <= nrstart
                startgrid[rr,cc]=1
            end
            cc = floor(Int, (xc-xmin)*ixrangece )+1
            rr = floor(Int, (yc-ymin)*iyrangere )+1
            if cc <= 0 || cc > ncend || rr <= 0 || rr > nrend
                break_j = j
                break
            elseif endgrid[rr,cc]==1
                if !(any(cc==tcc) && any(rr==trr))
                    break_j = j
                    break
                end
            else
                tcc=cc
                trr=rr
                endgrid[rr,cc]=1
            end
        end # end for j
        return v[1:break_j-1,:] 
    end  
     

    vertsout = Vector{Array{Float64,2}}()
    for idx in rc_list
        r, c = idx.I
        if !startgrid[r, c]
            startgrid[r, c] = 1
            x0 = xmin + (c - 0.5) * incstartx
            y0 = ymin + (r - 0.5) * incstarty

            # Correctly trace forward and backward
            vf = trace(x0, y0, itp_u, itp_v, stepsize, maxvert, xmin, xmax, ymin, ymax,  1)
            vb = trace(x0, y0, itp_u, itp_v, stepsize, maxvert, xmin, xmax, ymin, ymax, -1)

            vfo = process_trace(vf)
            vbo = process_trace(vb) 

            second_half = vbo[end:-1:1, :]
            first_half  = vfo[2:end, :]
            full_line = vcat(second_half, first_half)
            push!(vertsout, full_line)
            push!(vertsout, [NaN NaN])

        end # if startgrid
    end # end idx
    if isempty(vertsout)
        return reshape(Float64[], 0, 2)
    else
        return reduce(vcat, vertsout)
    end
    #return vcat(vertsout[1:end-1]...)
end 



"""
    get_streamlines(xx, yy, uu, vv;
                min_density::Real=1.0,
                max_density::Real=5.0)

Compute evenly-spaced streamlines on the grid `(xx,yy)` with vector field `(uu,vv)`.
If `compute_dist` is true, also return nearest-neighbor distances.
"""
function get_streamlines(xx, yy, uu, vv;
                     min_density::Real=1.0,
                     max_density::Real=5.0)

    
    if xx isa AbstractVector && yy isa AbstractVector
        xc, yc = xx, yy
        xg, yg = meshgrid(xc, yc)
    else
        xc = xx[1, :]; yc = yy[:, 1]
        xg, yg = xx, yy
    end
    @assert size(xg) == size(yg) == size(uu) == size(vv)
    @assert min_density > 0 && max_density > 0

    xy = get_stream_xy(xg, yg, uu, vv, xc, yc, min_density, max_density)
    return xy
end

