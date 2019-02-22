# Cameras for radiation pressure computation

"""
    ray_parallel(radius::Real)

Generate a ray along the -Z axis at a random point on a disk of given radius in the XY
plane.

"""
function ray_parallel(radius::Real)
    return ray_parallel(radius, rand(), rand())
end


"""
    ray_parallel(radius::Real, u1::Real, u2::Real)

Generate a ray along the -Z axis at a random point on a disk of given radius in the XY
plane.

This is the deterministic version of the function, taking two numbers `u1` and `u2` which
should be uniform random numbers in the [0, 1] interval to generate uniformly sampled
points on the disk.

"""
function ray_parallel(radius::Real, u1::Real, u2::Real)
    r = sqrt(u1) * radius
    phi = 2π * u2
    ray_origin = Point3(r * cos(phi), r * sin(phi), radius + 1.0)
    return Ray(ray_origin, -Z_AXIS)
end



"""
    ray_from_point(radius::Real, distance::Real, u1::Real, u2::Real)

Generate a Ray from a point which is the given distance on the Z axis, at a random point
on a sphere of given radius around the origin.

"""
function ray_from_point(radius::Real, distance::Real)
    return ray_from_point(radius, distance, rand(), rand())
end


"""
    ray_from_point(radius::Real, distance::Real, u1::Real, u2::Real)

Generate a Ray from a point which is the given distance on the Z axis, at a point on a
sphere of given radius around the origin.

This is the deterministic version of the function, taking two numbers `u1` and `u2` which
should be uniform random numbers in the [0, 1] interval to generate uniformly sampled
points in the cone from the distance to the visible cross section of the sphere.

"""
function ray_from_point(radius::Real, distance::Real, u1::Real, u2::Real)
    cos_th = sqrt(1.0 - (radius / distance)^2)
    z = cos_th + (1.0 - cos_th) * u1
    phi = 2π * u2
    r = sqrt(1.0 - z^2)
    direction = Vector3(r*cos(phi), r*sin(phi), -z)
    return Ray(Point3(0.0, 0.0, distance), direction)
end













