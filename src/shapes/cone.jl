
"""
    struct Cone <: Shape

Cone shape model.

"""
struct Cone <: Shape
    id::Int64
    h::Float64
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Cone() = Cone(1, 1.0, Transformation(), Transformation())
Cone(id::Integer) = Cone(id, 1.0, Transformation(), Transformation())
Cone(T::Transformation) = Cone(1, 1.0, T, inv(T))
Cone(id::Integer, T::Transformation) = Cone(id, 1.0, T, inv(T))
Cone(id::Integer, h::Real, T::Transformation) = Cone(id, h, T, inv(T))



function cone_between_points(id::Integer, P1::Point3, r1::Real, P2::Point3, r2::Real)
    if r2 ≈ r1
        error("Cone degenerate into cylinder (TODO: implement)")
    end
    reversed = r2 > r1
    if reversed
        r1, r2 = r2, r1
    end
        
    s = norm(P2 - P1)
    h = 1.0 - r2/r1
    zscale = s * r1 / (r1 - r2)
    
    flip = reversed ? translation(Z_AXIS * h) * rotation(X_AXIS, π) : Transformation()
    scale = scaling(r1, r1, zscale)
    rot = rotate_z_to(P2 - P1)
    trans = translation(P1...)
    
    C = Cone(id, h, trans*rot*scale*flip)
    
    return C
end




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
    v1 = Point3(1.0, 1.0, cone.h)
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
    Cone(P.id, P.h, T2, inv(T2))
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
    
    if (0.0 <= P1.z <= cone.h && tnear >= r.tmin)
        P = P1
        t = tnear
    elseif (0.0 <= P2.z <= cone.h && tfar <= r.tmax)
        P = P2
        t = tfar
    else
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
    
    # Compute hit points (on dual cone)
    P1 = ray(tnear)
    P2 = ray(tfar)

    # Return true if either point is on the "actual" cone surface
    return (0.0 <= P1.z <= cone.h && tnear > r.tmin) || (0.0 <= P2.z <= cone.h && tfar < r.tmax)
end
