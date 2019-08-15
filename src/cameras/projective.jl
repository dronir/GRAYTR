
"""
ProjectiveCamera can be 'regular' or orthographic depending on the projection.
"""
struct ProjectiveCamera{F<:Film} <: Camera
    film::F
    camera_to_world::Transformation
    raster_to_camera::Transformation
    camera_to_raster::Transformation
    camera_to_screen::Transformation
    screen_to_camera::Transformation
    dxcam::Vector3
    dycam::Vector3
    lensr::Float64
    focald::Float64
end

"""
    ProjectiveCamera(cam2world, projection, window, lensR, focal_depth, film)
    
Construct a ProjectiveCamera.
"""
function ProjectiveCamera(cam2w::Transformation, cam2screen::Transformation, window::Array{Float64,1}, lensr::Float64, focald::Float64, f::Film)
    screen2raster = (scaling(f.resX, f.resY, 1.0) 
                   * scaling(1.0/(window[2] - window[1]), 1.0/(window[4] - window[3]), 1.0) 
                   * translation(-window[1], -window[3], 0.0))
    raster2cam = inv(cam2screen) * inv(screen2raster)
    dxcam = raster2cam(Vector3(1,0,0))
    dycam = raster2cam(Vector3(0,1,0))
    ProjectiveCamera(f, cam2w, raster2cam, inv(raster2cam), cam2screen, inv(cam2screen), dxcam, dycam, lensr, focald)
end



"""
    OrthographicCamera(cam2world, window, lensR, focal_depth, film)

Construct an orthographic camera.
"""
function OrthographicCamera(cam2w::Transformation, window::Array{Float64,1}, lensr::Float64, focald::Float64, f::Film)
    ProjectiveCamera(cam2w, orthographic(0.0, 1.0), window, lensr, focald, f)
end


"""
    OrthographicCamera(cam2world, Xwidth, Ywidth, lensR, focal_depth, film)

Construct an orthographic camera.
"""
function OrthographicCamera(cam2w::Transformation, Xwidth::Real, Ywidth::Real, lensr::Float64, focald::Float64, f::Film)
    ProjectiveCamera(cam2w, orthographic(0.0, 1.0), [-Xwidth/2, Xwidth/2, -Ywidth/2, Ywidth/2], lensr, focald, f)
end


"""
    orthographic(znear::Real, zfar::Real)
    
Computes the orthographic projection transformation given near and far distances.
"""
function orthographic(znear::Real, zfar::Real)
    scaling(1.0, 1.0, 1.0/(zfar - znear)) * translation(0.0, 0.0, -znear)
end


"""
    generate_ray(C::ProjectiveCamera, sample::CameraSample)

Generate a ray leaving the camera based on the sample.
"""
function generate_ray(C::ProjectiveCamera, sample::CameraSample)
    pras = Point3(sample.imgX, sample.imgY, 0.0)
    pcam = C.raster_to_camera(pras)
    R = Ray(pcam, Vector3(0,0,1), 0.0, Inf, 1)
    return 1.0, C.camera_to_world(R)
end


"""
    count_tasks(camera::ProjectiveCamera, nCores::Integer)

A heuristic function to count the number of rendering tasks given the number of processor
cores and the resolution of the output image. This number is either four tasks per core,
or the number of 16x16 pixel blocks in the output, whichever is higher..

"""
function count_tasks(camera::ProjectiveCamera, nCores::Integer)
    npix = camera.film.resX * camera.film.resY
    nTasks = max(4*nCores, div(npix, 256))
    return round_pow2(nTasks)
end

