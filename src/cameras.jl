

struct ProjectiveCamera <: Camera
    film::Film
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
    screen2raster = (scaling(f.xres, f.yres, 1.0) 
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
    return C.camera_to_world(R)
end



struct ImageFilm{T<:Filter} <: Film
    resX::Int64
    resY::Int64
    filter::T
    pixels::Array{Float64,3}
    filtertable::Array{Float64,2}
end

const FILTERTABLE_SIZE = 16

function make_filtertable(f::Filter)
    ftbl = zeros(Float64, (FILTERTABLE_SIZE, FILTERTABLE_SIZE))
    for j = 1:FILTERTABLE_SIZE
        fy = (j + 0.5) * f.ywidth / FILTERTABLE_SIZE
        for i = 1:FILTERTABLE_SIZE
            fx = (i + 0.5) * f.xwidth / FILTERTABLE_SIZE
            ftbl[i,j] = evaluate(f, fx, fy)
        end
    end
    ftbl
end

ImageFilm(x::Integer, y::Integer, f::Filter) = ImageFilm(x, y, f, zeros(Float64, (x,y,3)), 
                                                         make_filtertable(f))

function add_sample(F::ImageFilm, sample::Sample, L::Spectrum)
    dimgX = sample.imgX - 0.5
    dimgY = sample.imgY - 0.5
    x0 = ceil(Int64, dimgX - F.filter.xwidth)
    x1 = floor(Int64, dimgX + F.filter.xwidth)
    y0 = ceil(Int64, dimgY - F.filter.ywidth)
    y1 = floor(Int64, dimgY + F.filter.ywidth)
    if (x1-x0) < 0 || (y1-y0) < 0
        return
    end
    xyz = spectrum_to_xyz(L)
    for i = x0:x1
        for j = y0:y0
            # TODO: unoptimized way, without using filter table
            w = evaluate(F.filter, (i-dimgX), (j-dimgY))
            F.pixels[i,j,:] += xyz
        end
    end
end
