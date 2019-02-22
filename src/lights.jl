


################################
# Visibility Tester is a ray from point to light source

# VisibilityTester between two points (with epsilons)
function VisibilityTester(p1::Point3, p2::Point3, eps1::Float64, eps2::Float64)
    distance = norm(p2-p1)
    return Ray(p1, (p2-p1)/distance, eps1, distance*(1.0 - eps2), 1)
end

# VisibilityTester from a point to infinitely far away
function VisibilityTester(origin::Point3, direction::Vector3, eps::Float64)
    return Ray(origin, direction, eps, Inf, 1)
end



################################
# Point light

struct PointLight{S<:Spectrum} <: LightSource
    position::Point3
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
end

# Minimal constructor
function PointLight(L::Spectrum, l2w::Transformation)
    PointLight(l2w(Point3(0)), L, l2w, inv(l2w))
end

isdelta(L::PointLight) = true
power(L::PointLight) = 4π * L.intensity

function sample_L(light::PointLight, p::Point3)
    distance = norm(light.position - p)
    L = light.intensity / distance^2
    wi = (light.position - p) / distance
    pdf = 1.0
    vis = VisibilityTester(p, light.position, 2e-5, 0.0)
    return L, wi, pdf, vis
end

direct(L::PointLight) = true
background(L::PointLight) = nolight

function generate_ray(light::PointLight, bounds::BoundingSphere, S::Sample)
    ray = ray_from_point(bounds.radius, S.imgX, S.imgY)
    return light.light_to_world(ray)
end


################################
# Distanct light source

struct DistantLight{S<:Spectrum} <: LightSource
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
end

function DistantLight(L::Spectrum, l2w::Transformation)
    DistantLight(L, l2w, inv(l2w))
end

function sample_L(light::DistantLight, p::Point3)
    return light.intensity, NEG_Z_AXIS, 1.0, VisibilityTester(p, NEG_Z_AXIS, 2e-5)
end

direct(L::DistantLight) = true
background(L::DistantLight) = nolight


function generate_ray(light::DistantLight, bounds::BoundingSphere, S::Sample)
    ray = ray_parallel(bounds.radius, S.imgX, S.imgY)
    return light.light_to_world(ray)
end




################################
# Background light source

struct Background{S<:Spectrum} <: LightSource
    intensity::S
end

sample_L(light::Background, p::Point3) = nolight

direct(L::Background) = false
background(L::Background) = L.intensity


# TODO: Area light
#
struct AreaLight <: LightSource end


################################
# Disk Light (TODO)

struct DiskLight{S<:Spectrum} <: LightSource 
    direction::Vector3
    radius::Float64
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
end


function DiskLight(direction::Vector3, radius::Float64, intensity::Spectrum)
    rot_z = rotate_z_to(direction)
    return DiskLight(direction, radius, intensity, rot_z, inv(rot_z))
end

background(L::DiskLight) = nolight
direct(L::DiskLight) = true


sample_L(light::DiskLight, p::Point3) = sample_L(light, p, rand(), rand())

function sample_L(light::DiskLight, p::Point3, u1::Real, u2::Real)
    transform = light.light_to_world
    direction = transform(disk_light_sample(light.radius, u1, u2))
    return light.intensity, light.direction, 1.0, VisibilityTester(p, direction, 2e-5)
end

function disk_light_sample(radius::Real, u1::Real, u2::Real)
    cos_th = cos(radius)
    z = cos_th + (1.0 - cos_th) * u1
    phi = 2π * u2
    r = sqrt(1.0 - z^2)
    return Vector3(r*cos(phi), r*sin(phi), z)
end


