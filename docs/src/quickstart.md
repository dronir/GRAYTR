
# Quick start


```
function full_computation()
    
    # Create spectrum with flat white color
    whitespec = SampledSpectrum(300, 800, [1.0 for i = 1:256])
    
    # Create spectrum with reddish colour
    redspec = SampledSpectrum(300, 800, collect(range(0.0, 1.0, length=256)))

    # Create a Lambertian surface material with the white spectrum
    mat = Lambert(whitespec)
    
    # Create sphere shape located at origin
    Tsph = translation(0.0, 0.0, 0.0)
    sph = Sphere(1, 1.0, Tsph)
    
    # Create open cylinder
    Tcyl = translation(2.0, 2.0, 0.0) * scaling(1.0, 1.0, 2.0)
    cyl = Cylinder(Tcyl)
    
    # Create primitives combining a shape and a material
    sph_primitive = GeometricPrimitive(sph, mat, nothing, 1)
    cyl_primitive = GeometricPrimitive(cyl, mat, nothing, 2)
    
    # Create list of primitives
    primitives = GeometricPrimitive[sph_primitive, cyl_primitive]
    
    # Generate bounding box hierarchy from primitives
    stuff = BVHAccelerator(primitives)
    
    # Create a distant light source with the red spectrum
    p_light = DistantLight(redspec, IDENTITY_TRANSFORM)
    
    # Create scene from the bounding box hierarchy and list of light sources
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
    @time render(scene, image_cam, whitted, sampler ; debug=true)
    
    # Write the output
    write_image(F, "test.png")

    return true
end

```



