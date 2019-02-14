

function full_pressure()
    
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
    lights = LightSource[p_light]
    
    # Create scene with bounding box hierarchy and light source
    scene = Scene(stuff, lights)
    
    # Set up pixel sampler with 9 rays per pixel (3 rays in both X and Y directions)
    res = 256
    sampler = StratifiedSampler(res, res, 3)
    
    # Declare an integrator with maximum ray depth of 1, using 10000 rays per light source
    # and initialize empty arrays to store the force and torque from each light source.
    Nrays = 10000
    forces = [Vector3(0) for i in lights]
    torques = [Vector3(0) for i in lights]
    integrator = PressureIntegrator(Nrays, 1, forces, torques)
    
    # Run the renderer
    @time render(scene, integrator, sampler)
    

    return true
    
    
end

@testset "Full pressure computation" begin
    @test full_pressure()
end


