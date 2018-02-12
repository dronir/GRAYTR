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
include("accelerators.jl")
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
mat = MatteMaterial(RGBSpectrum(1.0, 1.0, 1.0))

T1 = rotation(Y_AXIS, π/5) * rotation(X_AXIS, pi/6)
sph1 = Cylinder(T1)
sphP1 = GeometricPrimitive(sph1, mat, Nullable{AreaLight}(), 1)

T2 = translation(0.0, 0.0, 1.0)
sph2 = Disk(T1 * T2)
sphP2 = GeometricPrimitive(sph2, mat, Nullable{AreaLight}(), 1)

sph3 = Disk(T1 * rotation(X_AXIS, 0.0))
sphP3 = GeometricPrimitive(sph3, mat, Nullable{AreaLight}(), 1)

primitives = Primitive[sphP1, sphP2, sphP3]

println("Generating bounding box hierarchy...")
stuff = BVHAccelerator(primitives)

    

function make_light(col::Array, pos::Transformation)
    light_col = RGBSpectrum(col...) * 400
    return PointLight(light_col, pos, 1)
end


light1 = make_light([1.0, 0.0, 0.0], rotation(Z_AXIS, 0.0) * translation(10, 10, -10.0))
light2 = make_light([0.0, 1.0, 0.0], rotation(Z_AXIS, 4pi/5) * translation(10, 10, -10.0))
light3 = make_light([0.0, 0.0, 1.0], rotation(Z_AXIS, 6pi/5) * translation(10, 10, -10.0))
bg = Background(RGBSpectrum(10.0, 10.0, 10.0))

scene = Scene(stuff, LightSource[light1, light2, light3, bg])




# Make a camera
resX = 512
resY = 512
F = ImageFilm(resX, resY, TriangleFilter(1, 1))
shift = translation(0.0, 0.0, -4.0)
Cam = OrthographicCamera(shift, [-2.0, 2, -2, 2], 0.0, 0.0, F)

# Make the renderer
sampler = StratifiedSampler(1, resX, 1, resY, 2, 2, true)
whitted = WhittedIntegrator(1)
renderer = SamplerRenderer(sampler, Cam, whitted, DummyVolumeIntegrator())

# Run the renderer
println("Starting render")
tic()
render(renderer, scene)
toc()

write_txt(F)

#Profile.print()

end # module
