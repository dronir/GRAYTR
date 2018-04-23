

struct SamplerRenderer{S<:Sampler, C<:Camera, I<:SurfaceIntegrator, V<:VolumeIntegrator} <: Renderer
    sampler::S
    camera::C
    surf_integrator::I
    vol_integrator::V
end


struct SamplerRendererTask{S<:Scene, R<:Renderer}
    scene::S
    renderer::R
#    sampler::Sampler
    number::Int64
    count::Int64
end

# TODO: this is a dummy verion
enqueue_and_run{S<:SamplerRendererTask}(Tasks::Array{S,1}) = [run(task) for task in Tasks]



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
    
    println("Running $nTasks render tasks on $(nprocs()) cores...")
    enqueue_and_run(tasks)
end

# Round an integer up to nearest power of two
round_pow2(n::Integer) = 2^convert(typeof(n), ceil(log(2, n)))

# Heuristically compute number of rendering tasks
function count_tasks(nCores::Integer, resX::Integer, resY::Integer)
    npix = resX * resY
    nTasks = max(4*nCores, div(npix, 256))
    return round_pow2(nTasks)
end


# The core functionality really starts here!
# This runs a SamplerRendererTask, generating a ray, computing the intensity along that ray,
# and adding it to the camera film.
function run(task::SamplerRendererTask)
    maybe_sampler = get_subsampler(task.renderer.sampler, task.number, task.count)
    isnull(maybe_sampler) && return nothing
    subsampler = get(maybe_sampler)
    
    max_samples = 0 # TODO
    state = 0
    while true
        samples, state = get_samples(subsampler, state)
        length(samples) == 0 && break
        intensities = Array{Spectrum}(size(samples))

        # generate camera rays
        for i = 1:length(samples)
            # find camera ray for sample[i]
            weight, ray = generate_ray(task.renderer.camera, samples[i])
            # evaluate radiance along ray
            if weight > 0.0
                Li, maybe_isect = intensity(task.renderer, task.scene, ray, samples[i])
                Ls = weight*Li
            else
                Li = NoLight()
                maybe_isect = Nullable{Intersection}()
            end
            sample = samples[i]
            if !isblack(Ls)
                if uses_isect(task.renderer.camera.film)
                    add_sample!(task.renderer.camera.film, sample, Ls, maybe_isect)
                else
                    add_sample!(task.renderer.camera.film, sample, Ls)
                end
            end
        end
    end
    return nothing
end

function intensity(renderer::SamplerRenderer, scene::Scene, r::Ray, sample::Sample)
    maybe_isect = intersect(r, scene)
    if !isnull(maybe_isect)
        # Ray hits a scene object. Get its contribution from surface integrator.
        isect = get(maybe_isect)
        Li = intensity(renderer.surf_integrator, renderer, scene, isect, r, sample)
    else
        # Ray doesn't hit any scenery. Add contribution from background light sources.
        Li = sum(background(light) for light in scene.lights)
    end
    # TODO: add volume integrator contribution
    return Li, maybe_isect
end

