module PrettyStreamlines


using Interpolations: interpolate, Gridded, Linear
using Random: shuffle!

include("Streamlines.jl")
include("PlotsRecipe.jl")

export compute_streamlines


end # module PrettyStreamlines
