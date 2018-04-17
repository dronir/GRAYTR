include("geometry.jl")

module Raytracer

using Geometry

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


function main()
    # Make something to look at
    mat = MatteMaterial(Lambert(RGBSpectrum(1.0, 1.0, 1.0)))
    
    T1 = rotation(Y_AXIS, π/5) * rotation(X_AXIS, pi/6)
    cyl = Cylinder(T1)
    cylP = GeometricPrimitive(cyl, mat, Nullable{AreaLight}(), 1)
    
    T2 = translation(0.0, 0.0, 1.0)
    cap1 = Disk(T1 * T2)
    capP1 = GeometricPrimitive(cap1, mat, Nullable{AreaLight}(), 1)
    
    cap2 = Disk(T1 * rotation(X_AXIS, pi))
    capP2 = GeometricPrimitive(cap2, mat, Nullable{AreaLight}(), 1)
    
    Tsph = translation(0.0, 1.0, 2.5)
    sph = Sphere(1, 1.5, Tsph)
    sphP = GeometricPrimitive(sph, mat, Nullable{AreaLight}(), 1)
    
    primitives = Primitive[cylP, capP1, capP2, sphP]
    
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
    
    p_light = DistantLight(-Z_AXIS, RGBSpectrum(1,1,1)*1, rotation(X_AXIS, pi/6), 1)
    
    
    scene = Scene(stuff, LightSource[light2, p_light])
    
    
    
    
    # Make a camera
    resX = 512
    resY = 512
    F = ImageFilm(resX, resY, TriangleFilter(1, 1))
    shift = translation(0.0, 0.0, -2.0)
    
    DF = DelayFilm(532.0, resX, resY, 1000, 1.0, 5.0)
    
    window = [-2.0, 2.5, -2.0, 2.5]
    image_cam = OrthographicCamera(shift, window, 0.0, 0.0, F)
    delay_cam = OrthographicCamera(shift, window, 0.0, 0.0, DF)
    
    
    # Make the renderer
    sampler = StratifiedSampler(1, resX, 1, resY, 3, 3, true)
    whitted = WhittedIntegrator(1)
    image_renderer = SamplerRenderer(sampler, image_cam, whitted, DummyVolumeIntegrator())
    delay_renderer = SamplerRenderer(sampler, delay_cam, whitted, DummyVolumeIntegrator())
    
    
    # Run the renderer
    println("Starting image render")
    @time render(image_renderer, scene)
    #println("Starting delay render")
    #@time render(delay_renderer, scene)
    
    write_image(F)
    #write_txt(DF, "test.txt")
    
    #Profile.print()
end # function main()

#main()

end # module
