
# Dummy volume integrator for testing
struct DummyVolumeIntegrator <: VolumeIntegrator end
preprocess(D::DummyVolumeIntegrator) = true


struct WhittedIntegrator <: SurfaceIntegrator
    maxdepth::Int64
end

preprocess(W::WhittedIntegrator) = true

# This inner part of loop separated so that it can be dispatched based on the type of 'light':
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
        return (brdf_value .* light_spectrum) .* abs(dot(w0, dg.n))
    end
end

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


