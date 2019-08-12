
export LightSource, DistantLight, PointLight, DiskLight


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
# No background for direct light sources

"""
    background(L::DirectLight)

Give the background contribution of the given light source. Return `nolight` for a
`DirectLight`.

"""
@inline background(L::DirectLight) = nolight




################################
# Point light


"""
    PointLight{S<:Spectrum} <: DirectLight

Point light source. It's located at a given point and radiates isotropically. The intensity
value given is interpreted as the intensity in Watts/m^2 at a distance of 1 meter.

# Fields

- `position::Point3`, the location of the light source.
- `intensity::S`, the intensity at a distance of 1 meter in Watts/m^2.
- `light_to_world::Transformation`, transformation to shift the location of the light
   source.
- `world_to_light::Transformation`, inverse transformation of the above.

"""
struct PointLight{S<:Spectrum} <: DirectLight
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
    pdf = 1.0
    vis = VisibilityTester(p, light.position, 2e-5, 0.0)
    return L, pdf, vis
end


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
    DistantLight{S<:Spectrum} <: DirectLight
    
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
struct DistantLight{S<:Spectrum} <: DirectLight
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
- A weight coefficient which is always 1.0
- A ray from the light source to the given point.

"""
function sample_L(light::DistantLight, p::Point3)
    return light.intensity, 1.0, VisibilityTester(p, light.light_to_world(NEG_Z_AXIS), 1e-9)
end


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
    Background{S<:Spectrum} <: IndirectLight
    
Background light. If one of these is in the scene, every ray that escapes to infinity hits
this light source.

"""
struct Background{S<:Spectrum} <: IndirectLight
    intensity::S
end


"""
    sample_L(light::Background, p::Point3)

Generate a sample from light source. Return `nolight` because background lights do not
contribute through this process (they are not direct lights).

"""
sample_L(light::Background, p::Point3) = nolight


"""
    background(L::Background)
    
Returns the background light produced by the lightsource. Returns the intensity of a
`Background` source.

"""
background(L::Background) = L.intensity




###### ########################
# TODO: Area light
struct AreaLight <: LightSource end




################################
# Disk Light (TODO)

"""
    DiskLight{S<:Spectrum} <: DirectLight 

A light source which is a uniformly illuminated disk far away. Expressed in terms of its
angular radius. By default it's in the positive Z direction, and can be rotated with
a transformation.

# Fields:

- `radius::Float64`, angular radius (in radians) of the source
- `intensity::S`, 
- `light_to_world::Transformation`, transformation to adjust the direction of the source.
   Chosen to rotate the +Z axis into the desired direction of the light.
- `world_to_light::Transformation`, inverse transformation of the above

"""
struct DiskLight{S<:Spectrum} <: DirectLight 
    radius::Float64
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
end


"""
    DiskLight(direction::Vector3, radius::Float64, intensity::Spectrum)

Create new `DiskLight` in the desired direction, angular radius and intensity.

"""
function DiskLight(direction::Vector3, radius::Float64, intensity::Spectrum)
    rot_z = rotate_z_to(direction)
    return DiskLight(radius, intensity, rot_z, inv(rot_z))
end


"""
    background(L::DiskLight)
    
Returns the background light produced by the lightsource. Returns `nolight` for a
`DiskLight` source.

"""
background(L::DiskLight) = nolight


sample_L(light::DiskLight, p::Point3) = sample_L(light, p, rand(), rand())

function sample_L(light::DiskLight, p::Point3, u1::Real, u2::Real)
    direction = light.light_to_world(disk_light_sample(light.radius, u1, u2))
    return light.intensity, 1.0, VisibilityTester(p, direction, 2e-9)
end


"""
    disk_light_sample(radius::Real, u1::Real, u2::Real)

Create a vector which points from origin to a point on the surface of a disk with given
angular radius centered on the +Z axis. It is parametrized by two numbers, `u1` and `u2`,
which range from zero to one. Using uniform random numbers will result in the surface of
the disk to be uniformly sampled.

"""
function disk_light_sample(radius::Real, u1::Real, u2::Real)
    cos_th = cos(radius)
    z = cos_th + (1.0 - cos_th) * u1
    phi = 2π * u2
    r = sqrt(1.0 - z^2)
    return Vector3(r*cos(phi), r*sin(phi), z)
end


