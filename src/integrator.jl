
# Dummy volume integrator for testing
struct DummyVolumeIntegrator <: VolumeIntegrator end
preprocess(D::DummyVolumeIntegrator) = true


struct WhittedIntegrator <: SurfaceIntegrator
    maxdepth::Int64
end

preprocess(W::WhittedIntegrator) = true

# This inner part of loop separated so that it can be dispatched based on the type of 'light':
function inner_int(light::LightSource, p::Point3, bsdf::BSDF, wo::Vector3, n::Normal3, scene::Scene)
    !direct(light) && return nolight
    Li, wi, pdf, vis = sample_L(light, p)
    if isblack(Li) || pdf ≈ 0.0
        return nolight
    end
    if !intersectP(vis, scene)
        refl = evaluate(bsdf, wi, wo)
        return (refl .* Li) .* abs(dot(wi, n))
    end
    return nolight
end

function intensity(intgr::WhittedIntegrator, scene::Scene, isect::Intersection, 
                   ray::Ray, sample::Sample)
    # evaluate BSDF at intersect point
#    L = RGBSpectrum(0,0,0)
    bsdf = get_BSDF(isect)::BSDF
    p = bsdf.dgs.p
    n = bsdf.dgs.n
    wo = -ray.direction

    # compute emitted light if hit object is an area light
#    L += emission(iSect, w0)
    
    # add contribution of each light source
    L = nolight
    for light in scene.lights
        L += inner_int(light, p, bsdf, wo, n, scene)::Spectrum
    end
    return L
    
    # TODO: recursive for depth
#    if ray.depth < maxDepth
        # trace rays for transmission and reflection
#    end
end


