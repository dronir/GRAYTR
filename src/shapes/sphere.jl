
# Sphere type and related methods

"""
    Sphere <: Shape

Sphere shape model.

"""
struct Sphere <: Shape
    id::Int64
    radius::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Sphere() = Sphere(1, 1.0, false, Transformation(), Transformation())
Sphere(T::Transformation) = Sphere(1, 1.0, false, T, inv(T))
Sphere(r::Real, T::Transformation) = Sphere(1, r, false, T, inv(T))
Sphere(id::Int64, r::Real, T::Transformation) = Sphere(id, r, false, T, inv(T))
Sphere(id::Int64, r::Real, inverted::Bool, T::Transformation) = Sphere(id, r, inverted, T, inv(T))


"""
    can_intersect(s::Sphere)

Is the sphere an intersectable object? Always returns true.

"""
can_intersect(s::Sphere) = true

"""
    area(s::Sphere)

Compute surface area of a sphere (doesn't into account scaling...)

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
    Sphere(S.id, S.radius, S.inverted, T2, inv(T2))
end



"""
    shape_intersect(r::Ray, sph::Sphere)

Compute intersection geometry between given ray and sphere.

"""
function shape_intersect(r::Ray, sph::Sphere)
    ray = sph.world_to_obj(r)
    A = dot(ray.direction, ray.direction)
    B = 2.0 * dot(ray.direction, ray.origin)
    C = dot(ray.origin, ray.origin) - sph.radius^2
    hit, tnear, tfar = quadratic(A, B, C)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return nothing, NaN, NaN
    end
    t = tnear > r.tmin ? tnear : tfar
    P = ray(t)
    
    n = normalize(Normal3(P.x, P.y, P.z))
    s = normalize(Vector3(-2π*P.y, 2π*P.x, 0.0))
    
    DG = DifferentialGeometry(sph.obj_to_world(P), sph.obj_to_world(n), sph.obj_to_world(s))
    return DG, t, 5e-4 * t
end



"""
    intersectP(r::Ray, sph::Sphere)

Quick true/false intersection test between ray and sphere.

"""
function intersectP(r::Ray, sph::Sphere)
    ray = sph.world_to_obj(r)
    A = dot(ray.direction, ray.direction)
    B = 2.0 * dot(ray.direction, ray.origin)
    C = dot(ray.origin, ray.origin) - sph.radius^2
    hit, tnear, tfar = quadratic(A, B, C)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return false
    else
        return true
    end
end
