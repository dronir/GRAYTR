
"""
    struct Cone <: Shape

Cone shape model.

"""
struct Cone <: Shape
    id::Int64
    h::Float64
    hmax::Float64
    radius::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Cone() = Cone(1, 1.0, 1.0, 1.0, false, Transformation(), Transformation())
Cone(T::Transformation) = Cone(1, 1.0, 1.0, 1.0, false, T, inv(T))

"""
    can_intersect(s::Cone) = true

Is a paraboloid intersectable?

TODO: should this function be replaced by an abstract type like `IntersectableShape`?

"""
can_intersect(s::Cone) = true


"""
    obj_bounds(P::Cone)

The bounding box or a paraboloid in object space.

"""
function obj_bounds(P::Cone)
    v0 = Point3(-P.radius, -P.radius, P.h0)
    v1 = Point3(P.radius, P.radius, P.h1)
    return BoundingBox(v0, v1)
end


"""
    world_bounds(P::Cone)

The bounding box or a paraboloid transformed to world space.

"""
function world_bounds(P::Cone)
    P.obj_to_world(obj_bounds(P))
end


"""
    (T::Transformation)(P::Cone)

Apply a `Transformation` to a `Cone` object, return a new object.

"""
function (T::Transformation)(P::Cone)
    T2 = T * P.obj_to_world
    Cone(P.id, P.h0, P.h1, P.radius, P.inverted, T2, inv(T2))
end



function cone_coefs(ray::Ray, cone::Cone)
    A = cone.h^2 * (ray.direction.x^2 + ray.direction.y^2)  - cone.radius^2 * ray.origin.z^2
    B = 2*cone.h^2 * (ray.origin.x * ray.direction.x + ray.origin.y * ray.direction.y) - 2 * cone.radius^2 * ray.direction.z * (cone.h - ray.origin.z)
    C = cone.h^2 * (ray.origin.x^2 + ray.origin.y^2) - cone.radius^2 * (cone.h^2 - 2*cone.h*ray.origin.z + ray.origin.z^2)
    return A, B, C
end


function shape_intersect(r::Ray, cone::Cone)
    ray = cone.world_to_obj(r)
    
    hit, tnear, tfar = quadratic(cone_coefs(ray, cone)...)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return nothing, NaN, NaN
    end
    
    # Compute hit point, check if it's on our cone
    t = tnear > r.tmin ? tnear : tfar
    P = ray(t)
    if P.z > cone.h || P.z < 0.0 || P.z > cone.hmax
        return nothing, NaN, NaN
    end
    
    phi = atan(P.y, P.x)
    theta = atan(cone.radius, cone.h)
    
    n = Normal3(cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta))
    s = normalize(Vector3(cone.radius, 0.0, cone.h))
    
    n = flip_normal(n, ray.direction)
    
    DG = DifferentialGeometry(cone.obj_to_world(P), cone.obj_to_world(n), cone.obj_to_world(s))
    
    return DG, t, 5e-4 * t

end



"""
    intersectP(r::Ray, cone::Cone)

"Quick" true/false intersect test between given ray and paraboloid.

"""
function intersectP(r::Ray, cone::Cone)
    ray = cone.world_to_obj(r)
    
    hit, tnear, tfar = quadratic(cone_coefs(ray, cone)...)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return false
    end
    
    # Compute hit point, check if it's on our paraboloid.
    t = tnear > r.tmin ? tnear : tfar
    P = ray(t)
    if P.z > cone.h || P.z < 0.0 || P.z > cone.hmax
        return false
    else
        return true
    end
end