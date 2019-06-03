
"""
    struct Paraboloid <: Shape

Paraboloid shape model.

"""
struct Paraboloid <: Shape
    id::Int64
    h0::Float64
    h1::Float64
    radius::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Paraboloid() = Paraboloid(1, 0.0, 1.0, 1.0, false, Transformation(), Transformation())
Paraboloid(T::Transformation) = Paraboloid(1, 0.0, 1.0, 1.0, false, T, inv(T))

"""
    can_intersect(s::Paraboloid) = true

Is a paraboloid intersectable?

TODO: should this function be replaced by an abstract type like `IntersectableShape`?

"""
can_intersect(s::Paraboloid) = true


"""
    obj_bounds(P::Paraboloid)

The bounding box or a paraboloid in object space.

"""
function obj_bounds(P::Paraboloid)
    v0 = Point3(-P.radius, -P.radius, P.h0)
    v1 = Point3(P.radius, P.radius, P.h1)
    return BoundingBox(v0, v1)
end


"""
    world_bounds(P::Paraboloid)

The bounding box or a paraboloid transformed to world space.

"""
function world_bounds(P::Paraboloid)
    P.obj_to_world(obj_bounds(P))
end


"""
    (T::Transformation)(P::Paraboloid)

Apply a `Transformation` to a `Paraboloid` object, return a new object.

"""
function (T::Transformation)(P::Paraboloid)
    T2 = T * P.obj_to_world
    Paraboloid(P.id, P.h0, P.h1, P.radius, P.inverted, T2, inv(T2))
end



function shape_intersect(r::Ray, Par::Paraboloid)
    ray = Par.world_to_obj(r)
    A = Par.h1 * (ray.direction.x^2 + ray.direction.y^2)
    B = 2*Par.h1 * (ray.origin.x * ray.direction.x + ray.origin.y * ray.direction.y) - Par.radius^2 * ray.direction.z
    C = Par.h1 * (ray.origin.x^2 + ray.origin.y^2) - Par.radius^2 * ray.origin.z
    
    hit, tnear, tfar = quadratic(A, B, C)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return nothing, NaN, NaN
    end
    
    # Compute hit point, check if it's on our paraboloid.
    t = tnear > r.tmin ? tnear : tfar
    P = ray(t)
    if P.z > Par.h1 || P.z < Par.h0
        return nothing, NaN, NaN
    end
    
    r = sqrt(P.x^2 + P.y^2)
    
    n = normalize(Normal3(2r, 0.0, -1.0))
    s = normalize(Vector3(2r, 0.0, 1.0))
    
    n = flip_normal(n, ray.direction)
    
    DG = DifferentialGeometry(Par.obj_to_world(P), Par.obj_to_world(n), Par.obj_to_world(s))
    
    return DG, t, 5e-4 * t

end



"""
    intersectP(r::Ray, Par::Paraboloid)

"Quick" true/false intersect test between given ray and paraboloid.

"""
function intersectP(r::Ray, Par::Paraboloid)
    ray = Par.world_to_obj(r)
    A = Par.h1 * (ray.direction.x^2 + ray.direction.y^2)
    B = 2.0 * Par.h1 *(ray.origin.x * ray.direction.x + ray.origin.y * ray.direction.y) - Par.radius^2 * ray.direction.z
    C = Par.h1 * (ray.origin.x^2 + ray.origin.y^2) - Par.radius^2 * ray.origin.z
    
    hit, tnear, tfar = quadratic(A, B, C)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return false
    end
    
    # Compute hit point, check if it's on our paraboloid.
    t = tnear > r.tmin ? tnear : tfar
    P = ray(t)
    if P.z > Par.h1 || P.z < Par.h0
        return false
    else
        return true
    end
end