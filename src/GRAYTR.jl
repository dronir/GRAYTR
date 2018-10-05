
module GRAYTR

include("geometry.jl")
using .Geometry

using FileIO
using Images

include("abstracts.jl")
include("utility.jl")
include("rays.jl")
include("diffgeom.jl")
include("shape.jl")
include("samplers.jl")
include("bounding.jl")

include("spectrum.jl")
include("materials.jl")
include("scene.jl")
include("lights.jl")
include("primitives.jl")
include("accelerators.jl")
include("integrator.jl")
include("renderer.jl")
include("filters.jl")
include("cameras.jl")

include("parser.jl")

include("films/imagefilm.jl")
include("films/delayfilm.jl")


end #module
