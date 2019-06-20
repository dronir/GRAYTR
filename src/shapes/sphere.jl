
# Sphere type and related methods

"""
    Sphere <: Shape

Sphere shape model.

"""
struct Sphere <: Shape
    id::Int64
    radius::Float64
    angle::Float64
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Sphere() = Sphere(1, 1.0, π, Transformation(), Transformation())
Sphere(T::Transformation) = Sphere(1, 1.0, π, T, inv(T))
Sphere(r::Real, T::Transformation) = Sphere(1, r, π, T, inv(T))
Sphere(r::Real, angle::Real, T::Transformation) = Sphere(1, r, angle, T, inv(T))
Sphere(id::Int64, r::Real, T::Transformation) = Sphere(id, r, π, T, inv(T))
Sphere(id::Int64, r::Real, angle::Real, T::Transformation) = Sphere(id, r, angle, T, inv(T))


"""
    can_intersect(s::Sphere)

Is the sphere an intersectable object? Always returns true.

"""
can_intersect(s::Sphere) = true

"""
    area(s::Sphere)

Compute surface area of a sphere (doesn't into account scaling...)

TODO: fix for angle (pen and paper required)

"""
area(s::Sphere) = 4π * s.radius^2


"""
    obj_bounds(s::Sphere)

Return bounding box of the sphere in object coordinates.

"""
function obj_bounds(s::Sphere)
    v = Point3(s.radius, s.radius, s.radius)
    return BoundingBox(-v, v)
end


"""
    world_bounds(s::Sphere)

Return bounding box of the sphere in world coordinates.

"""
function world_bounds(s::Sphere)
    s.obj_to_world(obj_bounds(s))
end


"""
    (T::Transformation)(S::Sphere)

Apply a transformation to a sphere.

"""
function (T::Transformation)(S::Sphere)
    T2 = T * S.obj_to_world
    Sphere(S.id, S.radius, S.angle, T2, inv(T2))
end



function sphere_coefs(ray::Ray, sph::Sphere)
    A = dot(ray.direction, ray.direction)
    B = 2.0 * dot(ray.direction, ray.origin)
    C = dot(ray.origin, ray.origin) - sph.radius^2
    return A, B, C
end


"""
    shape_intersect(r::Ray, sph::Sphere)

Compute intersection geometry between given ray and sphere.

"""
function shape_intersect(r::Ray, sph::Sphere)
    ray = sph.world_to_obj(r)
    
    hit, tnear, tfar = quadratic(sphere_coefs(ray, sph)...)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return nothing, NaN, NaN
    end
    P1 = ray(tnear)
    P2 = ray(tfar)
    if P1.z / sph.radius < cos(sph.angle)
        if P2.z / sph.radius < cos(sph.angle)
            return nothing, NaN, NaN
        else
            P = P2
            t = tfar
        end
    else
        P = P1
        t = tnear
    end
    
    n = normalize(Normal3(P.x, P.y, P.z))
    s = normalize(Vector3(-2π*P.y, 2π*P.x, 0.0))
    
    n = flip_normal(n, ray.direction)
    
    DG = DifferentialGeometry(sph.obj_to_world(P), sph.obj_to_world(n), sph.obj_to_world(s))
    return DG, t, 5e-4 * t
end



"""
    intersectP(r::Ray, sph::Sphere)

Quick true/false intersection test between ray and sphere.

"""
function intersectP(r::Ray, sph::Sphere)
    ray = sph.world_to_obj(r)
   
    hit, tnear, tfar = quadratic(sphere_coefs(ray, sph)...)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return false
    end
    
    P1 = ray(tnear)
    P2 = ray(tfar)
    if P1.z / sph.radius < cos(sph.angle)
        if P2.z / sph.radius < cos(sph.angle)
            return false
        else
            return true
        end
    end
    
    return true
end
