
absract type Renderer end

struct SamplerRenderer <: Renderer
    sampler::Sampler
    camera::Camera
    surf_integrator::SurfaceIntegrator
    vol_integrator::VolumeIntegrator
end

abstract type Task end

struct SamplerRendererTask <: Task
    scene::Scene
    renderer::Renderer
    sample::Sample
    number::Int64
    count::Int64
end



# Render a given scene with a SamplerRenderer
function render(renderer::SamplerRenderer, scene::Scene)
    # integrator preprocessing
    preprocess(renderer.surf_integerator)
    preprocess(renderer.vol_integrator)
    # initialize sample storage
    # i.e. make an array of Samples to be filled?
    
    # create and launch tasks for rendering
    
    nTasks = count_tasks(NCORES, camera.film.resX, camera.film.resY)
    tasks = [SamplerRendererTask(scene, renderer, sample, n, nTasks) for n = 1:nTasks]
    enqueue_and_run(tasks)
    
    # cleanup and produce image
    write_image(renderer.camera.film)
end

# Round an integer up to nearest power of two
round_pow2(n::Integer) = 2^convert(typeof(n), ceil(log(2, n)))

# Heuristically compute number of rendering tasks
function count_tasks(nCores::Integer, resX::Integer, resY::Integer)
    npix = resX * resY
    nTasks = max(32*nCores, div(npix, 256))
    return round_pow2(nTasks)
end


function run(task::SamplerRendererTask)
    maybe_sampler = get_subsampler(task.sampler, task.number, task.count)
    isnull(maybe_sampler) && return nothing # What does it return?
    subsampler = get(maybe_sampler)
    
    max_samples = 0 # TODO
    array{Intersection}(max_samples)
    
    while true
        samples = get_samples(subsampler)
        size(samples) == 0 && break

        # generate camera rays
        for i = 1:size(samples)
            # find camera ray for sample[i]
            weight, ray = generate_raydifferential(task.renderer.camera)
            scaledifferentials(ray, 1/sqrt(sampler.samplesperpixel))
            # evaluate radiance along ray
            if weight ≈ 0
                Ls = weight * intensity(task.renderer, task.scene, task.ray, samples[i])
        end
        # TODO; report to Sampler
        # contribute to image
        add_sample(task.renderer.camera.film)
        
    end
end

function intensity(renderer::SamplerRenderer, r::RayDifferential, sample::Sample)
    Li = 0.0
    if intersect(r, renderer.scene)
        Li += intensity(renderer.surf_integerator, renderer, r, sample)
    else
        Li += sum(background(light) for light in renderer.scene.lights)
    end
    # TODO: add volume integrator contribution
    return Li
end


