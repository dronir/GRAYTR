

struct SamplerRenderer <: Renderer
    sampler::Sampler
    camera::Camera
    surf_integrator::SurfaceIntegrator
    vol_integrator::VolumeIntegrator
end


struct SamplerRendererTask
    scene::Scene
    renderer::Renderer
#    sampler::Sampler
    number::Int64
    count::Int64
end

# TODO: this is a dummy verion
enqueue_and_run(Tasks::Array{SamplerRendererTask,1}) = [run(task) for task in Tasks]



# Render a given scene with a SamplerRenderer
function render(renderer::SamplerRenderer, scene::Scene)
    # integrator preprocessing
    preprocess(renderer.surf_integrator)
    preprocess(renderer.vol_integrator)
    # initialize sample storage
    # i.e. make an array of Samples to be filled?
    
    # create and launch tasks for rendering
    camera = renderer.camera
    nTasks = count_tasks(nprocs(), camera.film.resX, camera.film.resY)
    tasks = [SamplerRendererTask(scene, renderer, n, nTasks) for n = 1:nTasks]
    enqueue_and_run(tasks)
    
    # cleanup and produce image
    write_image(renderer.camera.film)
end

# Round an integer up to nearest power of two
round_pow2(n::Integer) = 2^convert(typeof(n), ceil(log(2, n)))

# Heuristically compute number of rendering tasks
function count_tasks(nCores::Integer, resX::Integer, resY::Integer)
    npix = resX * resY
    nTasks = max(4*nCores, div(npix, 256))
    return round_pow2(nTasks)
end


function run(task::SamplerRendererTask)
    maybe_sampler = get_subsampler(task.renderer.sampler, task.number, task.count)
    isnull(maybe_sampler) && return nothing
    subsampler = get(maybe_sampler)
#    @printf("(%d, %d, %d, %d)\n", subsampler.xstart, subsampler.xend, subsampler.ystart, subsampler.yend)
    
    max_samples = 0 # TODO
    #Array{Intersection}(max_samples)
    state = 0
    while true
        samples, state = get_samples(subsampler, state)
        length(samples) == 0 && break
        Ls = RGBSpectrum(0.0, 0.0, 0.0)

        # generate camera rays
        for i = length(samples)
            # find camera ray for sample[i]
            weight, ray = generate_ray(task.renderer.camera, samples[i])
#            ray = scaledifferentials(ray, 1/sqrt(sampler.samplesperpixel))
            # evaluate radiance along ray
            if !(weight ≈ 0)
                Ls = weight * intensity(task.renderer, task.scene, ray, samples[i])
            end
        end
        # TODO; report to Sampler
        # contribute to image
        for sample in samples
            add_sample(task.renderer.camera.film, sample, Ls) # TODO
        end
        
    end
    return nothing
end

function intensity(renderer::SamplerRenderer, scene::Scene, r::Ray, sample::Sample)
    Li = RGBSpectrum(0,0,0)
    maybe_isect = intersect(r, scene)
    if !isnull(maybe_isect)
        isect = get(maybe_isect)
        Li += intensity(renderer.surf_integrator, renderer, scene, isect, r, sample)
#    else
        #Li += sum(background(light) for light in scene.lights)
    end
    # TODO: add volume integrator contribution
    return Li
end

