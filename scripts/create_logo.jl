using Plots, PrettyStreamlines


# — 2) set up a grid
x = y = -3:0.01:3

# Define the four logo dots: (x, y, color)
dots = [
    (0.0, 1, RGB(0.22, 0.596, 0.149)),   # green 
    (cosd(30), -sind(30), RGB(0.584, 0.345, 0.698)),   # purple 
    (cosd(150), -sind(150), RGB(0.796, 0.235, 0.2)),   # green  
]

x_coords  = [d[1] for d in dots]
y_coords  = [d[2] for d in dots]
col_colors = [d[3] for d in dots]


strengths = [1,1,1]
centers = [(d[1],d[2]) for d in dots]
ε = 0.1  

p = plot(xlims=(-1.5,1.5), 
        ylims=(-1.5,1.5),
        aspect_ratio = :equal,
        legend        = false,
        framestyle    = :none)



# Plot
scatter(
    x_coords, y_coords,
    markersize    = 72,
    markerstrokewidth = 0,
    markercolor   = col_colors,
    aspect_ratio  = :equal,
    xlims         = (-2, 2),
    ylims         = (-2, 2),
    grid          = false,
    framestyle    = :none,
    legend        = false,
)


function u(x,y)
    sum( -Γ*(y - yc)/((x - xc)^2 + (y - yc)^2 + ε^2) )
end

function v(x,y)
    sum(  Γ*(x - xc)/((x - xc)^2 + (y - yc)^2 + ε^2) )
end

U = 1 
a = 1 

function  up(x,y) 
    r = hypot(x,y-1)
    θ = atan(y-,x)
    if r < a
        return NaN
    else
        return U*(1-a^2 * cos(2θ)/r^2)
    end 
end 

function  vp(x,y) 
    r = hypot(x,y-1)
    θ = atan(y-1,x)
    if r < a
        return NaN
    else
        return -U*a^2 * sin(2θ)/r^2
    end
end 

streamlines(x, y, up, vp,aspect_ratio=:equal)