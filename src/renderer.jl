
abstract type RenderTask end


"""
    SamplerRendererTask{S<:Scene, C<:Camera, I<:SurfaceIntegrator}

A wrapper for the scene, camera and integrator as well as a task number and total count of
tasks. This information represents a piece of the full render job and is used to render
the correct pixels of the output.

"""
struct SamplerRendererTask{S<:Scene, C<:Camera, I<:SurfaceIntegrator, T<:Sampler} <: RenderTask
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

The computation is multi-threaded on a single processor. The number of threads is one by
default, and can be changed by setting the `JULIA_NUM_THREADS` environmental variable. A
lock is used to prevent multiple data from accessing the camera film simultaneously. In
many films, a ray can affect multiple pixels or bins, including ones "belonging" to a
different thread.

"""
function enqueue_and_run(Tasks::Array{S,1}) where S<:RenderTask
    lock = Threads.SpinLock()
    Threads.@threads for task in Tasks
        run(task, lock)
    end
end






function enqueue_debug(Tasks::Array{S,1}) where S<:RenderTask
    lock = Threads.SpinLock()
    for (t, task) in enumerate(Tasks)
        run(task, lock)
    end
end



"""
    render(scene::Scene, camera::Camera, integrator::SurfaceIntegrator, sampler::Sampler)

Render a given scene, as seen by a given camera, using given SurfaceIntegrator and Sampler.
"""
function render(scene::Scene, camera::Camera, integrator::SurfaceIntegrator, 
                sampler::Sampler ; debug=false)
    # initialize sample storage
    # i.e. make an array of Samples to be filled?
    
    # create and launch tasks for rendering
    nTasks = count_tasks(camera, nprocs())
    tasks = [SamplerRendererTask(scene, camera, integrator, sampler, n, nTasks) for n = 1:nTasks]
    
    if debug
        println("Running $(nTasks) render tasks in single-thread debug mode...")
        enqueue_debug(tasks)
    else
        println("Running $(nTasks) render tasks on $(Threads.nthreads()) threads...")
        enqueue_and_run(tasks)
    end
end



"""
    run(task::SamplerRendererTask, write_lock::Threads.AbstractLock)

Run a SamplerRendererTask. This is the core function of the raytracer: here we generate a
ray from the camera, trace its possible intersection with the scene, compute reflection of
all light sources at that point and add up the intensity to the camera film.

The thread lock given by the second argument, `write_lock` is used to lock the camera film
so it can only be accessed by one thread at a time. If a significant amount of time is
spent in the `add_sample!` function, multi-threading could lead to a performance drop. This
depends on the film.

"""
function run(task::SamplerRendererTask, write_lock::Threads.AbstractLock)
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
                    Ls = weight * intensity(task.integrator, task.scene, maybe_isect, 
                                            ray, samples[i])
                    lock(write_lock)
                    add_sample!(task.camera.film, samples[i], Ls, maybe_isect)
                    unlock(write_lock)
                end
            end
        end
    end
    return nothing
end






"""
    PressureRendererTask{S<:Scene, T<:Sampler} <: RenderTask

A structure that contains the information needed to render one part of the scene.

"""
struct PressureRendererTask{S<:Scene, T<:Sampler} <: RenderTask
    scene::S
    integrator::PressureIntegrator
    sampler::T
    nlight::Int64
    number::Int64
    count::Int64
    sphere::BoundingSphere
end


"""
    function render(scene::Scene, integrator::PressureIntegrator, sampler::Sampler)

Compute the radiation pressure and torque on a given scene, using given PressureIntegrator
and Sampler.

"""
function render(scene::Scene, integrator::PressureIntegrator, sampler::Sampler ; debug=false)
    # create and launch tasks for rendering
   
    N_lights = length(scene.lights)
    nbins = sampler.xend * sampler.yend
    rays_per_bin = sampler.xs * sampler.ys
    tasks_per_light = round_pow2(max(div(nbins * rays_per_bin, 1024), 4*nprocs()))
    N_tasks = N_lights * tasks_per_light
    tasks = PressureRendererTask[]
    sphere = BoundingSphere(scene.bounds)
    for (l, light) in enumerate(scene.lights)
        for n = 1:tasks_per_light
            push!(tasks, PressureRendererTask(scene, integrator, sampler, l,
                    n, tasks_per_light, sphere))
        end
        println()
    end
    
    if debug
        println("Running $(N_tasks) render tasks in single-thread debug mode...")
        enqueue_debug(tasks)
    else
        println("Running $(N_tasks) render tasks on $(Threads.nthreads()) threads...")
        enqueue_and_run(tasks)
    end
    
end





function run(task::PressureRendererTask, write_lock::Threads.AbstractLock)
    subsampler = get_subsampler(task.sampler, task.number, task.count)
    subsampler == nothing && return nothing
    Nsamples = subsampler.xs * subsampler.ys
    samples = Array{CameraSample}(undef, Nsamples)

    light = task.scene.lights[task.nlight]
    
    state = 0

    while !finished(subsampler, state)
        state = get_samples!(subsampler, state, samples)
        
        force = Vector3(0)
        torque = Vector3(0)
        
        for i = 1:length(samples)
            sample = normalize(samples[i], subsampler.xnorm, subsampler.ynorm)
            ray = generate_ray(light, task.sphere, sample)
            isect = intersect(ray, task.scene)
            if isect != nothing
                f, t = compute_pressure(task.integrator, light, isect, ray)
                force += f
                torque += t
            end
           task.integrator.counts[task.nlight] += 1
        end
        lock(write_lock)
        task.integrator.force[task.nlight] += force * cross_section(task.sphere)
        task.integrator.torque[task.nlight] += torque * cross_section(task.sphere)
        unlock(write_lock)
    end
end

