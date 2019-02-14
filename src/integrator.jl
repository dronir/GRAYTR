
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
    # evaluate BSDF at intersect point
    world_to_local = local_transformation(isect.geometry)

    # compute emitted light if hit object is an area light
#    L += emission(iSect, w0)
    
    # add contribution of each light source
    L = nolight
    for light in scene.lights
        L += inner_int(light, isect.geometry, isect.material, -ray.direction, world_to_local, scene)
    end
    return L
    
    # TODO: recursive for depth
#    if ray.depth < maxDepth
        # trace rays for transmission and reflection
#    end
end

"""
    inner_int(light::LightSource, dg::DifferentialGeometry, 
                 mat::BxDF, w1::Vector3, T::Transformation, scene::Scene)
                 
The inner loop of the `intensity` function, separated for better optimization.
"""
function inner_int(light::LightSource, dg::DifferentialGeometry, mat::BxDF, w1::Vector3, T::Transformation, scene::Scene)
    !direct(light) && return nolight
    light_spectrum, w0, pdf, light_ray = sample_L(light, dg.p)
    if isblack(light_spectrum) || pdf ≈ 0.0
        return nolight
    end
    if intersectP(light_ray, scene)
        return nolight
    else
        brdf_value = evaluate(mat, T, w1, w0)
        return (brdf_value * light_spectrum) * abs(dot(w0, dg.n))
    end
end





"""
    PressureIntegrator

A Whitted integrator that computes radiation pressure force and torque.
"""
struct PressureIntegrator <: SurfaceIntegrator
    n_rays::Int64
    maxdepth::Int64
    force::Array{Vector3,1}
    torque::Array{Vector3,1}
end


"""
    PressureIntegrator(n_rays::Integer)

Constructor for single-scattering integrator. Equivalent to `PressureIntegrator(n_rays, 1)`.

"""
PressureIntegrator(n_rays::Integer) = PressureIntegrator(n_rays, 1)



const TARGET_RAYS_PER_TASK = 1024

"""
    count_tasks(integrator::PressureIntegrator, scene::Scene, nprocs::Integer)

Return the number of PressureRendererTasks per light source for a given
`PressureIntegrator`, and number of processors `nprocs`.

Uses a heuristic which tries to keep the number of rays per task close to a constant 1024,
but also uses at least two tasks per processor core.

"""
function count_tasks(integrator::PressureIntegrator, nprocs::Integer)
    rays_per_light = integrator.n_rays
    m1 = div(rays_per_light, TARGET_RAYS_PER_TASK)
    m2 = 2*nprocs
    
    return max(m1, m2)
end





"""

"""
function compute_pressure(intgr::PressureIntegrator, scene::Scene, isect::Intersection, 
                   ray::Ray, sample::Sample)
    #
    

#    p = compute_pressure(isect.material, L, direction)
    f = Vector3(0)
    t = cross(f, isect.geometry.p)
    
    return Vector3(0), Vector3(0)
end
