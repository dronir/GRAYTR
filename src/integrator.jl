
# Dummy volume integrator for testing
struct DummyVolumeIntegrator <: VolumeIntegrator end

preprocess(D::DummyVolumeIntegrator) = true

struct WhittedIntegrator <: SurfaceIntegrator
    maxdepth::Int64
end

preprocess(W::WhittedIntegrator) = true

+(x::Void, y::Spectrum) = y

function inner_int(light::LightSource, p::Point3, bsdf, wo::Vector3, n::Normal3, scene::Scene)
    if !direct(light)
        return NoLight()
    end
    Li, wi, pdf, vis = sample_L(light, p)
    if isblack(Li) || pdf ≈ 0.0
        return NoLight()
    end
    refl = evaluate(bsdf, wi, wo)
    if !isblack(refl) && unoccluded(vis, scene)
        return (refl * Li) * abs(dot(wi, n))
    end
    return NoLight()
end

function intensity(intgr::WhittedIntegrator, rend::Renderer, scene::Scene, isect::Intersection, 
                   ray::Ray, sample::Sample)
    # evaluate BSDF at intersect point
#    L = RGBSpectrum(0,0,0)
    bsdf = get_BSDF(isect)
    p = bsdf.dgs.p
    n = bsdf.dgs.n
    wo = -ray.direction

    # compute emitted light if hit object is an area light
#    L += emission(iSect, w0)
    
    # add contribution of each light source
    L = NoLight()
    for light in scene.lights
        L += inner_int(light, p, bsdf, wo, n, scene) 
    end
    return L
    
    # TODO: recursive for depth
#    if ray.depth < maxDepth
        # trace rays for transmission and reflection
#    end
end

