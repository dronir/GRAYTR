


struct WhittedIntegrator <: SurfaceIntegrator
    maxdepth::Int64
end

function intensity(intgr::WhittedIntegrator, rend::Renderer, scene::Scene, isect::Intersection, 
                   r::Ray, sample::Sample)
    # evaluate BSDF at intersect point
    L = 0.0
    bsdf = get_BSDF(isect, ray)
    p = bsdf.dgshading.p
    n = bsdf.dgshading.n
    wo = -r.direction

    # compute emitted light if hit object is an area light
    L += emission(iSect, w0)
    
    # add contribution of each light source
    for light in scene.lights
        Li, wi, pdf, vis = sample_light(light)
        if isblack(Li) || pdf ≈ 0.0
            continue
        end
        f = bsdf(wi, wo)
        if !isblack(f) && !occluded(vis)
            L += f .* Li * abs(dot(wi, n))
        end
    end
    return L
    
    # TODO: recursive for depth
#    if ray.depth < maxDepth
        # trace rays for transmission and reflection
#    end
end

