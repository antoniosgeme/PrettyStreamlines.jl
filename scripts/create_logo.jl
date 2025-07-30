using Plots, PrettyStreamlines

dots = [
    (0.0, 1, RGB(0.22, 0.596, 0.149)),   # green 
    (cosd(30), -sind(30), RGB(0.584, 0.345, 0.698)),   # purple 
    (cosd(150), -sind(150), RGB(0.796, 0.235, 0.2)),   # green  
]

xc  = [d[1] for d in dots]
yc  = [d[2] for d in dots]
colors = [d[3] for d in dots]

x = y = collect(-3:0.01:3)
X = [i for j in y, i in x]
Y = [j for j in y, i in x]
Z = @. X + im * Y

V_inf = 1.0                     # freestream speed
a     = 3/4                     # cylinder radi
Γ     = [0.0, -0.0, 0.0]        # vortex strengths

α   = -π/2                      # freestream angle
zc  = xc .+ im .* yc            # comlpex cylinder centers
w   = V_inf * exp(-im*α) .* ones(size(Z))        # start with uniform flow

for i in eachindex(zc)
    z0 = zc[i]
    inside = abs.(Z .- z0) .< a
    Z[inside] .= NaN
    w .= w .- V_inf * a^2 * exp(im*α) ./ (Z .- z0).^2
    w .= w .- im*Γ[i] ./ (2π*(Z .- z0))
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

U = real.(w)
V = -imag.(w)

p = plot(
    aspect_ratio  = :equal,
    xlims         = extrema(x),
    ylims         = extrema(y),
    grid          = false,
    framestyle    = :none,
    legend        = false,
    dpi           = 600
)

t = 0:0.01:2π
for i = 1:3
    plot!(3/4*cos.(t).+xc[i],3/4*sin.(t).+yc[i],fill=true,c=colors[i])
end

pressure(x,y,u,v) = 1 - (u^2 + v^2)
streamlines!(x,y,U,V,
            min_density = 1.3,
            max_density = 10,
            color_by=pressure,
            cmap=:jet  ,
            lw=2,
            glyphs=true,
            unbroken=true,
            arrow_every=60,
)

display(p)

savefig("./assets/logo.png")