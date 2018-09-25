

module Raytracer

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



function satellite()
    whitespec = SampledSpectrum(300, 800, [1.0, 1.0, 1.0])
    mat1 = AshkhminShirleySingle(whitespec, 10.0, 0.8)
    mat1 = Lambert(whitespec)
    
    
    wing_rotation = translation(0.0, 0.0, 0.5) * rotation(X_AXIS, 0.0) * translation(0.0, 0.0, -0.5)
    shift_up = translation(0.0, 0.0, 1.0)

    shift_down = translation(0.0, 0.0, -0.5)
    overall_rotation = rotation(X_AXIS, pi/0.24) * rotation(Z_AXIS, pi/12) * shift_down

    cap1 = GeometricPrimitive(Disk(overall_rotation * shift_up), mat1, nothing, 1)    
    cap2 = GeometricPrimitive(Disk(overall_rotation), mat1, nothing, 1)
    cyl = GeometricPrimitive(Cylinder(overall_rotation), mat1, nothing, 1)
    
    delta = translation(1e-5, 0.0, 0.0)
    
    p1 = Point3(1, 0, 0)
    p2 = Point3(2, 0, 0)
    p3 = Point3(1, 0, 1)
    p4 = Point3(2, 0, 1)
    wing_A1 = GeometricPrimitive(Triangle(p1, p4, p2, overall_rotation*wing_rotation), mat1, nothing, 1)
    wing_A2 = GeometricPrimitive(Triangle(p1, p3, p4, overall_rotation*wing_rotation*delta), mat1, nothing, 1)
    
    p1 = Point3(-1, 0, 0)
    p2 = Point3(-2, 0, 0)
    p3 = Point3(-1, 0, 1)
    p4 = Point3(-2, 0, 1)
    wing_B1 = GeometricPrimitive(Triangle(p1, p4, p2, overall_rotation*wing_rotation), mat1, nothing, 1)
    wing_B2 = GeometricPrimitive(Triangle(p1, p3, p4, overall_rotation*wing_rotation*delta), mat1, nothing, 1)
    
    return GeometricPrimitive[cap1, cap2, cyl, wing_A1, wing_A2, wing_B1, wing_B2]
end



function main()
    # Make something to look at
    matspec = SampledSpectrum(300, 800, [0.3, 0.3, 0.3])
    whitespec = SampledSpectrum(300, 800, [1.0, 1.0, 1.0])
    mat1 = AshkhminShirleySingle(whitespec, 50.0, 0.8)
    mat2 = Lambert(whitespec)
    
    Tsph = translation(0.0, 0.0, 0.0)
    sph = Sphere(1, 1.0, Tsph)
    sphP = GeometricPrimitive(sph, mat2, nothing, 1)
    
    Ttri = Transformation()
    tri1 = Triangle(Ttri)
    triP1 = GeometricPrimitive(tri1, mat2, nothing, 1)
    
    lambert_sphere = GeometricPrimitive(Sphere(1, 1.0, Tsph), mat2, nothing, 1)
    lambert_disk = GeometricPrimitive(Disk(rotation(Z_AXIS, pi/2)), mat2, nothing, 1)
    
    
#    primitives = GeometricPrimitive[cylP, capP1, capP2]
    primitives = satellite()
    
    println("Generating bounding box hierarchy...")
    stuff = BVHAccelerator(primitives)
    
        
    
    #function make_light(col::Array, pos::Transformation)
    #    light_col = SampledSpectrum(300, 800, col)
    #    return PointLight(light_col, pos, 1)
    #end
    #
    #white_light = SampledSpectrum(300, 800, [1.0, 1.0, 1.0]*2)
    #
    #light1 = make_light([1.0, 0.0, 0.0], rotation(Z_AXIS, 0.0) * translation(10, 10, -10.0))
    #light2 = make_light([0.0, 1.0, 0.0], rotation(Z_AXIS, 4pi/5) * translation(10, 10, -10.0))
    #light3 = make_light([0.0, 0.0, 1.0], rotation(Z_AXIS, 6pi/5) * translation(10, 10, -10.0))
    #bg = Background(white_light)
    
    p_light = DistantLight(-Z_AXIS, SingleLine(532.0, 1.0), rotation(X_AXIS, 0.0))
      
    scene = Scene(stuff, LightSource[p_light])
    
    
    
    # Make a camera
    resX = 512
    resY = 512
    width = 2.5
    window = [-width, width, -width, width]
    
    F = ImageFilm(resX, resY, 0.75, TriangleFilter(1, 1))
    shift = translation(0.0, 0.0, -2.0)
    
    
    DF = DelayFilm(532.0, window, resX, resY, 1000, 0.0, 4.0)
    
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
    println("Starting delay render")
    @time render(delay_renderer, scene)
    
    write_image(F)
    write_txt(DF, "test.txt")
    
    #Profile.print()
end # function main()

#main()

end # module
