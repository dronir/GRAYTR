

function full_pressure()
    
    # Create material with flat white spectrum
    whitespec = SampledSpectrum(300, 800, [1.0, 1.0, 1.0])
    mat = Lambert(whitespec)
    
    # Create sphere shape
    Tsph = translation(0.0, 0.0, 0.0)
    sph = Sphere(1, 1.0, Tsph)
    disk = Disk(translation(2, 0, 0))
        
    # Create primitive combining shape and material
    sph_primitive = GeometricPrimitive(sph, mat, nothing, 1)
    dsk_primitive = GeometricPrimitive(disk, mat, nothing, 1)
    
    
    # Create list of primitives
    primitives = GeometricPrimitive[sph_primitive, dsk_primitive]
    
    # Generate bounding box hierarchy from primitives
    stuff = BVHAccelerator(primitives)
    
    # Create a distant light source with a single-line 532 nm spectrum
    p_light = DistantLight(SingleLine(532.0, 1.0), rotation(X_AXIS, 0.0))
    p_light2 = DistantLight(SingleLine(532.0, 1.0), rotation(X_AXIS, pi))
    lights = LightSource[p_light]
    
    # Create scene with bounding box hierarchy and light source
    scene = Scene(stuff, lights)
    
    # Set up pixel sampler with 9 rays per pixel (3 rays in both X and Y directions)
    res = 256
    sampler = StratifiedSampler(res, res, 5)
    
    # Declare an integrator with maximum ray depth of 1, using 100000 rays per light source
    # and initialize empty arrays to store the force and torque from each light source.
    forces = [Vector3(0) for i in lights]
    torques = [Vector3(0) for i in lights]
    counts = [0 for i in lights]
    integrator = PressureIntegrator(1, forces, torques, counts)
    
    # Run the renderer
    @time render(scene, integrator, sampler ; debug=true)
    
    forces = forces ./ counts
    torques = torques ./ counts
    
    return forces, torques
end

@testset "Full pressure computation" begin
    f, t = full_pressure()
    @test isapprox(f[1].z, -13/9*π  - 5/3 * π ; atol=5)
end


