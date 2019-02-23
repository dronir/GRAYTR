


################################
# Visibility Tester is a ray from point to light source

"""
    VisibilityTester(p1::Point3, p2::Point3, eps1::Float64, eps2::Float64)

Returns a ray to test intersection on the line segment from `p1` to `p2` with the given ray
epsilon values.

"""
function VisibilityTester(p1::Point3, p2::Point3, eps1::Float64, eps2::Float64)
    distance = norm(p2-p1)
    return Ray(p1, (p2-p1)/distance, eps1, distance*(1.0 - eps2), 1)
end


"""
    VisibilityTester(origin::Point3, direction::Vector3, eps::Float64)

Returns a ray to test intersection on the line segment from `origin` infinitely far in the
given direction, with a given ray epsilon.

"""
function VisibilityTester(origin::Point3, direction::Vector3, eps::Float64)
    return Ray(origin, direction, eps, Inf, 1)
end




################################
# Point light


"""
    PointLight{S<:Spectrum} <: LightSource

Point light source. It's located at a given point and radiates isotropically. The intensity
value given is interpreted as the intensity in Watts/m^2 at a distance of 1 meter.

# Fields

- `position::Point3`, the location of the light source.
- `intensity::S`, the intensity at a distance of 1 meter in Watts/m^2.
- `light_to_world::Transformation`, transformation to shift the location of the light
   source.
- `world_to_light::Transformation`, inverse transformation of the above.

"""
struct PointLight{S<:Spectrum} <: LightSource
    position::Point3
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
end


"""
    PointLight(L::Spectrum, l2w::Transformation)


Create a new `PointLight` with the given spectrum, located at the origin but shifted by
the given `Transformation`.

"""
function PointLight(L::Spectrum, l2w::Transformation)
    PointLight(l2w(Point3(0)), L, l2w, inv(l2w))
end


"""
    sample_L(light::PointLight, p::Point3)

Get an intensity sample from the given `DistantLight` and a ray from the light towards the
given point. The sample is always the same since the light source is a uniform parallel
light.

Returns:

- Intensity `Spectrum`
- Direction `Vector3` of the light.
- A weight coefficient which is always 1.0
- A ray from the light source to the given point.

"""
function sample_L(light::PointLight, p::Point3)
    distance = norm(light.position - p)
    L = light.intensity / distance^2
    wi = (light.position - p) / distance
    pdf = 1.0
    vis = VisibilityTester(p, light.position, 2e-5, 0.0)
    return L, wi, pdf, vis
end


"""
    direct(L::PointLight)

Is the given light a direct light source? Returns true for a `PointLight`.

"""
direct(L::PointLight) = true


"""
    background(L::PointLight)

Give the background contribution of the given light source. Return `nolight` for a
`PointLight`.

"""
background(L::PointLight) = nolight


"""
    generate_ray(light::PointLight, bounds::BoundingSphere, S::Sample)

Returns a ray from the light source at the given `BoundingSphere`, parametrized by the
given random sample.

"""
function generate_ray(light::PointLight, bounds::BoundingSphere, S::Sample)
    ray = ray_from_point(bounds.radius, S.imgX, S.imgY)
    return light.light_to_world(ray)
end


################################
# Distanct light source

"""
    DistantLight{S<:Spectrum} <: LightSource
    
A point light source infinitely far away, i.e. parallel rays with a given intensity spectrum
interpreted as Watts / m^2.

By default the light source is "located" infinitely far away along the positive Z axis and
the light rays travel in the -Z direction.

# Fields

- `intensity::S`, the intensity of the light in Watts/m^2.
- `light_to_world::Transformation`, transformation to define the direction of the light
   source.
- `world_to_light::Transformation`, inverse transformation of the above.


"""
struct DistantLight{S<:Spectrum} <: LightSource
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
end


"""
    function DistantLight(L::Spectrum, l2w::Transformation)

Create a new `DistantLight` given the intensity spectrum (in units of Watts / m^2) and
a `Transformation` used to rotate the direction of the light.
    
"""
function DistantLight(L::Spectrum, l2w::Transformation)
    DistantLight(L, l2w, inv(l2w))
end


"""
    sample_L(light::DistantLight, p::Point3)

Get an intensity sample from the given `DistantLight` and a ray from the light towards the
given point. The sample is always the same since the light source is a uniform parallel
light.

Returns:

- Intensity `Spectrum`
- Direction `Vector3` of the light.
- A weight coefficient which is always 1.0
- A ray from the light source to the given point.

"""
function sample_L(light::DistantLight, p::Point3)
    return light.intensity, NEG_Z_AXIS, 1.0, VisibilityTester(p, NEG_Z_AXIS, 2e-5)
end


"""
    direct(L::DistantLight)

Is the given light a direct light? Returns `true` for a `DistantLight`.

"""
direct(L::DistantLight) = true


"""
    background(L::DistantLight)

Returns the background light produced by the lightsource. Returns `nolight` for a
`DistantLight`.

"""
background(L::DistantLight) = nolight


"""
    generate_ray(light::DistantLight, bounds::BoundingSphere, S::Sample)

Returns a ray from the light source at the given `BoundingSphere`, parametrized by the
given random sample.

"""
function generate_ray(light::DistantLight, bounds::BoundingSphere, S::Sample)
    ray = ray_parallel(bounds.radius, S.imgX, S.imgY)
    return light.light_to_world(ray)
end





################################
# Background light source

"""
    Background{S<:Spectrum} <: LightSource
    
Background light. If one of these is in the scene, every ray that escapes to infinity hits
this light source.

"""
struct Background{S<:Spectrum} <: LightSource
    intensity::S
end


"""
    sample_L(light::Background, p::Point3)

Generate a sample from light source. Return `nolight` because background lights do not
contribute through this process (they are not direct lights).

"""
sample_L(light::Background, p::Point3) = nolight


"""
    direct(L::Background)

Is this light source direct? Returns false for `Background` source.

"""
direct(L::Background) = false


"""
    background(L::Background)
    
Returns the background light produced by the lightsource. Returns the intensity of a
`Background` source.

"""
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


