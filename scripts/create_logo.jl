using Plots, PrettyStreamlines

# Define the four logo dots: (x, y, color)
dots = [
    (0.0, 1, RGB(0.22, 0.596, 0.149)),   # green 
    (cosd(30), -sind(30), RGB(0.584, 0.345, 0.698)),   # purple 
    (cosd(150), -sind(150), RGB(0.796, 0.235, 0.2)),   # green  
]

xc  = [d[1] for d in dots]
yc  = [d[2] for d in dots]
colors = [d[3] for d in dots]

# set up a grid
x = y = collect(-5:0.05:5)
X = [i for j in y, i in x]
Y = [j for j in y, i in x]
Z = @. X + im * Y


# Initialize plot

p = plot(
    aspect_ratio  = :equal,
    xlims         = (-2, 2),
    ylims         = (-2, 2),
    grid          = false,
    framestyle    = :none,
    legend        = false,
)

t = 0:0.01:2π
for i = 1:3
    plot!(3/4*cos.(t).+xc[i],3/4*sin.(t).+yc[i],fill=true,c=colors[i])
end

    # --- parameters ---
V_inf = 1.0                     # freestream speed
a     = 3/4                     # cylinder radius
Γ     = [0.0, -0.0, 0.0]        # vortex strengths (one per cylinder)

# --- precompute ---
α   = -π/2                      # freestream from +∞ downwards
zc  = xc .+ im .* yc            # cylinder centers
w   = V_inf * exp(-im*α) .* ones(size(Z))        # start with uniform flow

# --- loop over each cylinder i ---
for i in eachindex(zc)
    z0 = zc[i]

    # 1) mask out interior of cylinder i
    inside = abs.(Z .- z0) .< a
    Z[inside] .= NaN

    # 2) image‐doublet for uniform flow around cylinder i
    w .= w .- V_inf * a^2 * exp(im*α) ./ (Z .- z0).^2

    # 3) the real vortex at z0
    w .= w .- im*Γ[i] ./ (2π*(Z .- z0))

    # 4) images of the *other* vortices in cylinder i
    #= for j in eachindex(zc)
        if j == i
            continue
        end
        # location of vortex j
        z_vj    = zc[j]
        # inversion in cylinder i
        z_image = z0 + a^2/ conj(z_vj - z0)
        w .=  .+ im*Γ[j] ./ (2π*(Z .- z_image))
    end =#
end

# --- return Cartesian components ---
U = real.(w)
V = -imag.(w)

streamlines!(x,y,U,V)

display(p)