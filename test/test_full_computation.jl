
"""

    full_computation()

A test function which sets up an entire scene, with an object and a light source, as well
as an `ImageFilm` camera, renders the image and writes it to disk.

"""
function full_computation()
    
    # Create material with flat white spectrum
    whitespec = SampledSpectrum(300, 800, [1.0 for i = 1:256])
    redspec = SampledSpectrum(300, 800, collect(range(0.0, 1.0, length=256)))
    mat = Lambert(whitespec)
    
    # Create sphere shape
    Tsph = translation(0.0, 0.0, 0.0)
    sph = Sphere(1, 1.0, Tsph)
    
    # Create open cylinder
    Tcyl = translation(2.0, 1.0, 0.0) * scaling(1.0, 1.0, 2.0)
    cyl = Cylinder(Tcyl)
    
    # Create primitive combining shape and material
    sph_primitive = GeometricPrimitive(sph, mat, nothing, 1)
    cyl_primitive = GeometricPrimitive(cyl, mat, nothing, 2)
    
    # Create list of primitives
    primitives = GeometricPrimitive[sph_primitive, cyl_primitive]
    
    # Generate bounding box hierarchy from primitives
    stuff = BVHAccelerator(primitives)
    
    # Create a distant light source with a single-line 532 nm spectrum
    p_light = DistantLight(redspec, rotation(Z_AXIS, π/8) * rotation(X_AXIS, π/2))
    
    # Create scene with bounding box hierarchy and light source
    scene = Scene(stuff, LightSource[p_light])
    
    # Make a camera
    resX = 1024
    resY = 1024
    widthX = 5
    widthY = 5
    F = ImageFilm(resX, resY, 1.0, TriangleFilter(1, 1))
    camera_position = rotation(Y_AXIS, π/2) * translation(0.0, 0.0, -3.0)
    image_cam = OrthographicCamera(camera_position, widthX, widthY, 0.0, 0.0, F)
    
    # Set up pixel sampler with 9 rays per pixel (3 rays in both X and Y directions)
    sampler = StratifiedSampler(resX, resY, 3)
    
    # Declare an integrator with maximum ray depth of 1
    whitted = WhittedIntegrator(1)
    
    # Run the renderer
    @time render(scene, image_cam, whitted, sampler ; debug=true)
    
    # Write the output
    write_image(F, "test.png")

    return true
end


@testset "Full computation" begin
    @test full_computation()
    @test isfile("test.png")
end



