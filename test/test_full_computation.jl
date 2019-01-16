
"""

    full_computation()

A test function which sets up an entire scene, with an object and a light source, as well
as an `ImageFilm` camera, renders the image and writes it to disk.

"""
function full_computation()
    
    # Create material with flat white spectrum
    whitespec = SampledSpectrum(300, 800, [1.0, 1.0, 1.0])
    mat = Lambert(whitespec)
    
    # Create sphere shape
    Tsph = translation(0.0, 0.0, 0.0)
    sph = Sphere(1, 1.0, Tsph)
    
    # Create primitive combining shape and material
    sph_primitive = GeometricPrimitive(sph, mat, nothing, 1)
    
    # Create list of primitives
    primitives = GeometricPrimitive[sph_primitive]
    
    # Generate bounding box hierarchy from primitives
    stuff = BVHAccelerator(primitives)
    
    # Create a distant light source with a single-line 532 nm spectrum
    p_light = DistantLight(-Z_AXIS, SingleLine(532.0, 1.0), rotation(X_AXIS, 0.0))
    
    # Create scene with bounding box hierarchy and light source
    scene = Scene(stuff, LightSource[p_light])
    
    # Make a camera
    resX = 1024
    resY = 1024
    width = 2.5
    window = [-width, width, -width, width]
    F = ImageFilm(resX, resY, 0.75, TriangleFilter(1, 1))
    shift = translation(0.0, 0.0, -2.0)
    image_cam = OrthographicCamera(shift, window, 0.0, 0.0, F)
    
    # Set up pixel sampler with 9 rays per pixel (3 rays in both X and Y directions)
    sampler = StratifiedSampler(resX, resY, 3)
    
    # Declare an integrator with maximum ray depth of 1
    whitted = WhittedIntegrator(1)
    
    # Run the renderer
    @time render(scene, image_cam, whitted, sampler)
    
    # Write the output
    write_image(F, "test.png")

    return true
end


@testset "Full computation" begin
    @test full_computation()
    @test isfile("test.png")
end



