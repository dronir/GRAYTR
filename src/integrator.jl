
# Dummy volume integrator for testing
struct DummyVolumeIntegrator <: VolumeIntegrator end

preprocess(D::DummyVolumeIntegrator) = true

struct WhittedIntegrator <: SurfaceIntegrator
    maxdepth::Int64
end

preprocess(W::WhittedIntegrator) = true

function intensity(intgr::WhittedIntegrator, rend::Renderer, scene::Scene, isect::Intersection, 
                   ray::Ray, sample::Sample)
    # evaluate BSDF at intersect point
    L = RGBSpectrum(0,0,0)
    bsdf = get_BSDF(isect)
    p = bsdf.dgs.p
    n = bsdf.dgs.n
    wo = -ray.direction

    # compute emitted light if hit object is an area light
#    L += emission(iSect, w0)
    
    # add contribution of each light source
    for light in scene.lights
        if !direct(light)
            continue
        end
        Li, wi, pdf, vis = sample_L(light, p)
        if isblack(Li) || pdf ≈ 0.0
            continue
        end
        refl = evaluate(bsdf, wi, wo)
        if !isblack(refl) && unoccluded(vis, scene)
            L += refl .* Li * abs(dot(wi, n))
        end
    end
    return L
    
    # TODO: recursive for depth
#    if ray.depth < maxDepth
        # trace rays for transmission and reflection
#    end
end

