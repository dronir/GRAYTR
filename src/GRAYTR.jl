
module GRAYTR

include("geometry.jl")
using .Geometry

using LinearAlgebra
using Random
using Distributed
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

export SampledSpectrum, SingleLine
export GeometricPrimitive, Sphere, Disk, Cylinder, Triangle
export AshkhminShirleySingle, Lambert, LommelSeeliger
export Transformation, rotation, translation, scaling
export DistantLight, PointLight
export Scene, BVHAccelerator
export ImageFilm, DelayFilm, PhotometricFilm
export TriangleFilter, BoxFilter
export OrthographicCamera
export StratifiedSampler, WhittedIntegrator
export render
export write_image, write_txt

end #module
