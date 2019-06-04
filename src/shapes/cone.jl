
"""
    struct Cone <: Shape

Cone shape model.

"""
struct Cone <: Shape
    id::Int64
    top::Float64
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Cone() = Cone(1, 1.0, Transformation(), Transformation())
Cone(id::Integer) = Cone(id, 1.0, Transformation(), Transformation())
Cone(T::Transformation) = Cone(1, 1.0, T, inv(T))
Cone(id::Integer, T::Transformation) = Cone(id, 1.0, T, inv(T))


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
function obj_bounds(cone::Cone)
    v0 = Point3(-1.0, -1.0, 0.0)
    v1 = Point3(1.0, 1.0, cone.top)
    return BoundingBox(v0, v1)
end


"""
    world_bounds(P::Cone)

The bounding box or a paraboloid transformed to world space.

"""
function world_bounds(cone::Cone)
    cone.obj_to_world(obj_bounds(cone))
end


"""
    (T::Transformation)(P::Cone)

Apply a `Transformation` to a `Cone` object, return a new object.

"""
function (T::Transformation)(P::Cone)
    T2 = T * P.obj_to_world
    Cone(P.id, P.top, T2, inv(T2))
end



function cone_coefs(ray::Ray)
    A = ray.direction.x^2 + ray.direction.y^2  - ray.direction.z^2
    B = 2 * (ray.origin.x * ray.direction.x + ray.origin.y * ray.direction.y - ray.origin.z * ray.direction.z + ray.direction.z)
    C = (ray.origin.x^2 + ray.origin.y^2 - ray.origin.z^2 + 2*ray.origin.z - 1)
    return A, B, C
end


function shape_intersect(r::Ray, cone::Cone)
    ray = cone.world_to_obj(r)
    
    hit, tnear, tfar = quadratic(cone_coefs(ray)...)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return nothing, NaN, NaN
    end
    
    # Compute hit point, check if it's on our cone
    P1 = ray(tnear)
    P2 = ray(tfar)
    if P1.z < cone.top && P2.z < cone.top
        t = tnear > r.tmin ? tnear : tfar
        P = ray(t)
    elseif P1.z < cone.top
        t = tnear
        P = P1
    elseif P2.z < cone.top
        t = tfar
        P = P2
    else
        return nothing, NaN, NaN
    end
    if P.z > 1.0 || P.z > cone.top || P.z < 0.0
        return nothing, NaN, NaN
    end
    
    phi = atan(P.y, P.x)
    
    n = Normal3(cos(phi)*sqrt(2)/2, sin(phi)*sqrt(2)/2, sqrt(2)/2)
    s = normalize(Vector3(1.0-P.x, 1.0-P.y, 1.0-P.z))
    
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
    
    hit, tnear, tfar = quadratic(cone_coefs(ray)...)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return false
    end
    
    # Compute hit point, check if it's on our actual cone (not on the dual cone opening
    # upwards from z=1).
    P1 = ray(tnear)
    P2 = ray(tfar)
    
    # If both hit points are on the right cone, take the closer one, otherwise take
    # whichever is, or return false if neither is.
    if P1.z < cone.top && P2.z < cone.top
        t = tnear > r.tmin ? tnear : tfar
        P = ray(t)
    elseif P1.z < cone.top
        P = P1
    elseif P2.z < cone.top
        P = P2
    else
        return false
    end
    
    if P.z > 1.0 || P.z > cone.top || P.z < 0.0
        return false
    else
        return true
    end
end