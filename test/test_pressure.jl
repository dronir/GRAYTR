





function full_pressure()
    
    # Create material with flat white spectrum
    whitespec = SampledSpectrum(300, 800, [1.0, 1.0, 1.0])
    mat = Lambert(whitespec)
    
    # Create sphere shape
    sph = Sphere(1, 1.0, translation(0.0, 1.0, 0.0))
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



function setup_disk_pressure(angle::Real)
    whitespec = SampledSpectrum(300, 800, [1.0, 1.0, 1.0])
    mat = Lambert(whitespec)
    
    # Create sphere shape
    disk = Disk(translation(2, 0, 0))
        
    # Create primitive combining shape and material
    dsk_primitive = GeometricPrimitive(disk, mat, nothing, 1)
    
    
    # Create list of primitives
    primitives = GeometricPrimitive[dsk_primitive]
    
    # Generate bounding box hierarchy from primitives
    stuff = BVHAccelerator(primitives)
    
    # Create a distant light source with a single-line 532 nm spectrum
    # and intensity of 2.0 Watts / m^2 / sr
    rot = rotation(X_AXIS, angle)
    p_light = DiskLight(rot(Z_AXIS), 0.01, SingleLine(532.0, 2.0))
    lights = LightSource[p_light]
    
    # Create scene with bounding box hierarchy and light source
    scene = Scene(stuff, lights)
    
    # Set up pixel sampler with 9 rays per pixel (5 rays in both X and Y directions)
    res = 256
    sampler = StratifiedSampler(res, res, 5)
    
    # Declare an integrator with maximum ray depth of 1, using 100000 rays per light source
    # and initialize empty arrays to store the force and torque from each light source.
    forces = [Vector3(0) for i in lights]
    torques = [Vector3(0) for i in lights]
    counts = [0 for i in lights]
    integrator = PressureIntegrator(1, forces, torques, counts)
    
    return scene, integrator, sampler
end


@testset "Disk pressure with varying incidence" begin

    # Solid angle of the light source times intensity
    I0 = 2.0
    F0 = 2π * (1 - cos(0.01)) * I0

    # Test against analytical solution for various angles
    angles = [0, π/4, π/3, π/2]
    results = [-F0 * cos(angle) * (cos(angle) + 2//3) * π for angle in angles]
    
    for (angle, expected) in zip(angles, results)
        scene, integrator, sampler = setup_disk_pressure(angle)
        render(scene, integrator, sampler ; debug=true)
        force = integrator.force[1] / integrator.counts[1]
        
        @test isapprox(force.z, expected ; atol=1e-5)
    end

end





