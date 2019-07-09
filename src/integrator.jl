
"""
    WhittedIntegrator

The standard (and currently only) SurfaceIntegrator. Traces rays from surface point to light
sources and adds together their contributions. Recursion to ray depths greater than one is not
currently implemented yet.

"""
struct WhittedIntegrator <: SurfaceIntegrator
    "maxdepth, the maximum ray depth to trace recursively."
    maxdepth::Int64
end


"""
    intensity(integrator, scene, isect, ray, sample)

Return intensity given a WhittedIntegrator, a scene, intersection point, ray and camera
sample. 

# Arguments
- `integrator::WhittedIntegrator`: a WhittedIntegrator instance.
- `scene::Scene`: the Scene to trace illumination rays in.
- `isect::Intersection`: an Intersection, giving the local geometry at ray intersection.
- `ray::Ray`, the incident Ray (not used by this function but part of SurfaceIntegrator API)
- `sample::Sample`, a Sample (not used by this function but part of SurfaceIntegrator API)

"""
function intensity(intgr::WhittedIntegrator, scene::Scene, isect::Intersection, 
                   ray::Ray, sample::Sample)
    
    # add contribution of each light source
    L = nolight
    for light in scene.lights
        L = L .+ inner_int(light, isect.geometry, isect.material, -ray.direction, scene)
    end
    return L
    
    # TODO: recursive for depth
#    if ray.depth < maxDepth
        # trace rays for transmission and reflection
#    end
end


"""
    inner_int(light::DirectLight, dg::DifferentialGeometry, mat::BxDF, w1::Vector3, 
              T::Transformation, scene::Scene)
                 
The inner loop of the `intensity` function, separated for better optimization.

"""
function inner_int(light::DirectLight, dg::DifferentialGeometry, mat::BxDF, w1::Vector3, scene::Scene)
    light_spectrum, pdf, light_ray = sample_L(light, dg.p)
    if isblack(light_spectrum) || pdf ≈ 0.0
        return nolight
    end
    if intersectP(light_ray, scene)
        return nolight
    else
        w0 = -light_ray.direction
        brdf_value = evaluate(mat, dg.T, w1, w0)
        return (brdf_value .* light_spectrum) .* abs(dot(w0, dg.n))
    end
end


"""
    inner_int(light::IndirectLight, dg::DifferentialGeometry, 
              mat::BxDF, w1::Vector3, T::Transformation, scene::Scene)
                 
The inner loop of the `intensity` function, separated for better optimization. For an
`IndirectLight` this is always `nolight`.

"""
function inner_int(light::IndirectLight, dg::DifferentialGeometry, mat::BxDF, w1::Vector3, T::Transformation, scene::Scene)
    return nolight
end






"""
    PressureIntegrator

A Whitted integrator that computes radiation pressure force and torque.

"""
struct PressureIntegrator <: SurfaceIntegrator
    maxdepth::Int64
    force::Array{Vector3,1}
    torque::Array{Vector3,1}
    counts::Array{Int64,1}
end

"""

"""
function compute_pressure(intgr::PressureIntegrator, light::LightSource, isect::Intersection, 
                   ray::Ray)
    world_to_local = isect.geometry.T
    local_to_world = inv(world_to_local)
    incident_direction = world_to_local(-ray.direction)
    f = compute_pressure(isect.material, incident_direction, light.intensity)
    f = local_to_world(f)
    t = cross(f, isect.geometry.p)
    
    return f, t
end
