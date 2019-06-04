
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
include("brdf.jl")

include("ray_generators.jl")
include("lights.jl")
include("primitives.jl")
include("dumb_aggregate.jl")
include("accelerators.jl")
include("scene.jl")
include("integrator.jl")
include("renderer.jl")
include("filters.jl")

include("cameras/projective.jl")

include("parser.jl")

include("films/imagefilm.jl")
include("films/delayfilm.jl")


export Vector3, Point3, Normal3
export X_AXIS, Y_AXIS, Z_AXIS, NEG_X_AXIS, NEG_Y_AXIS, NEG_Z_AXIS, IDENTITY_TRANSFORM
export Transformation, rotation, translation, scaling

export SampledSpectrum, SingleLine
export GeometricPrimitive, Sphere, Disk, Cylinder, Triangle, Paraboloid, Cone
export AshkhminShirleySingle, Lambert, LommelSeeliger
export LightSource, DistantLight, PointLight
export Scene, BVHAccelerator
export ImageFilm, DelayFilm, PhotometricFilm
export TriangleFilter, BoxFilter
export OrthographicCamera
export StratifiedSampler, WhittedIntegrator, PressureIntegrator
export render
export write_image, write_txt, reset!

end #module
