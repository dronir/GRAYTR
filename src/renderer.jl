


"""
    SamplerRendererTask{S<:Scene, C<:Camera, I<:SurfaceIntegrator}

A wrapper for the scene, camera and integrator as well as a task number and total count of
tasks. This information represents a piece of the full render job and is used to render
the correct pixels of the output.

"""
struct SamplerRendererTask{S<:Scene, C<:Camera, I<:SurfaceIntegrator, T<:Sampler}
    scene::S
    camera::C
    integrator::I
    sampler::T
    number::Int64
    count::Int64
end



"""
    enqueue_and_run(Tasks::Array{S,1})

Add the given tasks to a queue and run them. Currently this is a dummy version which just
runs all the tasks one by one in a single thread.

TODO: This is where the code's parallelization should happen.
"""
function enqueue_and_run(Tasks::Array{S,1}) where S<:SamplerRendererTask
    for task in Tasks
        run(task)
    end
end



"""
    render(scene::Scene, camera::Camera, integrator::SurfaceIntegrator, sampler::Sampler)

Render a given scene, as seen by a given camera, using given SurfaceIntegrator and Sampler.
"""
function render(scene::Scene, camera::Camera, integrator::SurfaceIntegrator, 
                sampler::Sampler)
    # initialize sample storage
    # i.e. make an array of Samples to be filled?
    
    # create and launch tasks for rendering
    nTasks = count_tasks(nprocs(), camera.film.resX, camera.film.resY)
    tasks = [SamplerRendererTask(scene, camera, integrator, sampler, n, nTasks) for n = 1:nTasks]
    
    println("Running $nTasks render tasks on $(nprocs()) cores...")
    enqueue_and_run(tasks)
end


"""
    round_pow2(n::Integer)

Round an integer up to nearest power of two.
"""
round_pow2(n::Integer) = 2^convert(typeof(n), ceil(log(2, n)))


"""
    count_tasks(nCores::Integer, resX::Integer, resY::Integer)

A heuristic function to count the number of rendering tasks given the number of processor
cores and the resolution of the output image. This number is either four tasks per core,
or the number of 16x16 pixel blocks in the output, whichever is higher..
"""
function count_tasks(nCores::Integer, resX::Integer, resY::Integer)
    npix = resX * resY
    nTasks = max(4*nCores, div(npix, 256))
    return round_pow2(nTasks)
end


"""
    run(task::SamplerRendererTask)

Run a SamplerRendererTask. This is the core function of the raytracer: here we generate a
ray from the camera, trace its possible intersection with the scene, compute reflection of
all light sources at that point and add up the intensity to the camera film.
"""
function run(task::SamplerRendererTask)
    subsampler = get_subsampler(task.sampler, task.number, task.count)
    subsampler == nothing && return nothing
    
    max_samples = 0 # TODO
    state = 0
    
    Nsamples = subsampler.xs * subsampler.ys
    samples = Array{CameraSample}(undef, Nsamples)
    
    # TODO: This would more elegantly be an iterator, by defining iterate(subsampler, state)
    # and iterate(subsampler). But in that case the recycling of the `samples` array may
    # not work, which would lead to a performance drop. Experimentation needed.
    while true
        finished(subsampler, state) && break
        state = get_samples!(subsampler, state, samples)

        # generate camera rays
        for i = 1:length(samples)
            # find camera ray for sample[i]
            weight, ray = generate_ray(task.camera, samples[i])
            # evaluate radiance along ray and add it to camera film
            if weight > 0.0
                maybe_isect = intersect(ray, task.scene)
                if maybe_isect != nothing
                    Ls = weight * intensity(task.integrator, task.scene, maybe_isect, ray, samples[i])
                    if uses_isect(task.camera.film)
                        add_sample!(task.camera.film, samples[i], Ls, maybe_isect)
                    else
                        add_sample!(task.camera.film, samples[i], Ls)
                    end
                end
            end
        end
    end
    return nothing
end


