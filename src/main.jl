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

T1 = translation(0.1, 0.0, 0.0)
sph1 = Sphere(1, 1.0, T1)
sphP1 = GeometricPrimitive(sph1, mat, Nullable{AreaLight}(), 1)

T2 = translation(1.4, 1.0, -1.0)
sph2 = Sphere(1, 0.5, T2)
sphP2 = GeometricPrimitive(sph2, mat, Nullable{AreaLight}(), 1)

T3 = translation(-1.0, 0.0, -2.0)
sph3 = Sphere(1, 0.5, T3)
sphP3 = GeometricPrimitive(sph3, mat, Nullable{AreaLight}(), 1)

primitives = Primitive[sphP1, sphP2, sphP3]

println("Generating bounding box hierarchy...")
stuff = BVHAccelerator(primitives)


for node in stuff.nodes
    println(node)
end
    

function make_light(col::Vector3, pos::Point3)
    light_col = RGBSpectrum(cos) * 100
    light_pos = translation(pos)
    return PointLight(light_col, light_pos, 1)
    
end


light_col = RGBSpectrum(1.0, 1.0, 1.0) * 300
light_pos = translation(10.0, 10.0, -10.0)
light = PointLight(light_col, light_pos, 1)
bg = Background(RGBSpectrum(0.0, 0.0, 10.0))

scene = Scene(stuff, LightSource[light, bg])




# Make a camera
resX = 512
resY = 512
F = ImageFilm(resX, resY, TriangleFilter(0.5, 0.5))
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

write_bwtxt(F)

#Profile.print()

end # module
