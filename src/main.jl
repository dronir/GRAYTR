include("geometry.jl")

module Raytracer

using Geometry

using FileIO
using Images

include("abstracts.jl")
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
include("integrator.jl")
include("renderer.jl")
include("filters.jl")
include("cameras.jl")


# Make something to look at
T = Transformation()
sph = Sphere(1, 1.0, T)
mat = MatteMaterial(RGBSpectrum(1.0, 0.0, 0.0))
stuff = GeometricPrimitive(sph, mat, Nullable{AreaLight}(), 1)

light_col = RGBSpectrum(1.0, 1.0, 1.0) * 2
light_pos = translation(4.0, 4.0, -4.0)
light = PointLight(light_col, light_pos, 1)
bg = Background(RGBSpectrum(0.0, 0.0, 0.2))

scene = Scene(stuff, LightSource[light, bg])

# Make a camera
resX = 256
resY = 256
F = ImageFilm(resX, resY, TriangleFilter(2, 2))
shift = translation(0.0, 0.0, -2.0)
Cam = OrthographicCamera(shift, [-1.2, 1.2, -1.2, 1.2], 0.0, 0.0, F)

# Make the renderer
sampler = StratifiedSampler(1, resX, 1, resY, 4, 4, true)
whitted = WhittedIntegrator(1)
renderer = SamplerRenderer(sampler, Cam, whitted, DummyVolumeIntegrator())

# Run the renderer
println("Starting render")
tic()
render(renderer, scene)
toc()

#Profile.print()

end # module
