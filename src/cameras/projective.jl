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

function ProjectiveCamera(cam2w::Transformation, proj::Transformation, window::Array{Float64,1}, lensr::Float64, focald::Float64, f::Film)
    screen2raster = (scaling(f.resX, f.resY, 1.0) 
                   * scaling(1.0/(window[2] - window[1]), 1.0/(window[4] - window[3]), 1.0) 
                   * translation(-window[1], -window[3], 0.0))
    cam2screen = orthographic(0.0, 1.0)
    raster2screen = inv(screen2raster)
    raster2cam = inv(cam2screen) * raster2screen
    dxcam = raster2cam(Vector3(1,0,0))
    dycam = raster2cam(Vector3(0,1,0))
    ProjectiveCamera(f, cam2w, raster2cam, inv(raster2cam), cam2screen, inv(cam2screen), dxcam, dycam, lensr, focald)
end

function OrthographicCamera(cam2w::Transformation, window::Array{Float64,1}, lensr::Float64, focald::Float64, f::Film)
    ProjectiveCamera(cam2w, orthographic(0.0, 1.0), window, lensr, focald, f)
end

function orthographic(znear::Real, zfar::Real)
    scaling(1.0, 1.0, 1.0/(zfar - znear)) * translation(0.0, 0.0, -znear)
end

function generate_ray(C::ProjectiveCamera, sample::CameraSample)
    pras = Point3(sample.imgX, sample.imgY, 0.0)
    pcam = C.raster_to_camera(pras)
    R = Ray(pcam, Vector3(0,0,1), 0.0, Inf, 1)
    return 1.0, C.camera_to_world(R)
end
