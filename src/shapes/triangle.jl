
struct Triangle <: Shape
    id::Int64
    p1::Point3
    p2::Point3
    p3::Point3
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Triangle(T::Transformation) = Triangle(1, Point3(0,0,0), Point3(1,0,0), Point3(0,1,0), false, T, inv(T))


can_intersect(T::Triangle) = true
area(T::Triangle) = 0.5 * cross(T.p2 - T.p1, T.p3 - T.p1)

function obj_bounds(T::Triangle)
    return BoundingBox(T.p1, T.p2, T.p3)
end

function world_bounds(T::Triangle)
    T.obj_to_world(obj_bounds(T))
end

function shape_intersect(R::Ray, T::Triangle)
    ray = T.world_to_obj(R)
    
    ray = T.world_to_obj(R)
    
    e1 = T.p2 - T.p1
    e2 = T.p3 - T.p1
    
    s1 = cross(ray.direction, e2)
    divisor = dot(s1, e1)
    if divisor == 0.0
        return Nullable{DifferentialGeometry}(), NaN, NaN
    end
    invd = 1.0 / divisor
    
    d = ray.origin - T.p1
    b1 = dot(d, s1) * invd
    if b1 < 0.0 || b1 > 1.0
        return Nullable{DifferentialGeometry}(), NaN, NaN
    end
    
    s2 = cross(d, e1)
    b2 = dot(ray.direction, s2) * invd
    if b2 < 0.0 || b1+b2 > 1.0
        return Nullable{DifferentialGeometry}(), NaN, NaN
    end
    
    t = dot(e2, s2) * invd
    if t < ray.tmin || t > ray.tmax
        return Nullable{DifferentialGeometry}(), NaN, NaN
    end
    
    P = ray(t)
    
    u = b1
    v = b2
    
    du1 = 0.0
    du2 = 1.0
    dv1 = -1.0
    dv2 = -1.0
    dp1 = T.p1 - T.p3
    dp2 = T.p2 - T.p3
    determinant = du1*dv2 - dv1*du2
    invdet = 1.0 / determinant
    dpdu = ( dv2 * dp1 - dv1 * dp2) * invdet
    dpdv = (-du2 * dp1 + du1 * dp2) * invdet
    
    DG = Nullable(DifferentialGeometry(
        T.obj_to_world(P),
        u, v, T,
        T.obj_to_world(dpdu),
        T.obj_to_world(dpdv),
        Normal3(0,0,0),
        Normal3(0,0,0)
    ))
    return DG, t, 5e-4 * t
end


# Fast true/false intersection test between a ray and a triangle
function intersectP(R::Ray, T::Triangle)
    ray = T.world_to_obj(R)
    
    e1 = T.p2 - T.p1
    e2 = T.p3 - T.p1
    
    s1 = cross(ray.direction, e2)
    divisor = dot(s1, e1)
    if divisor == 0.0
        return false
    end
    invd = 1.0 / divisor
    
    d = ray.origin - T.p1
    b1 = dot(d, s1) * invd
    if b1 < 0.0 || b1 > 1.0
        return false
    end
    
    s2 = cross(d, e1)
    b2 = dot(ray.direction, s2) * invd
    if b2 < 0.0 || b1+b2 > 1.0
        return false
    end
    
    t = dot(e2, s2) * invd
    if t < ray.tmin || t > ray.tmax
        return false
    end
    
    return true
end




