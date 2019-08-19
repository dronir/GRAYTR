
struct Triangle <: Shape
    id::Int64
    p1::Point3
    p2::Point3
    p3::Point3
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Triangle(T::Transformation) = Triangle(1, Point3(0,0,0), Point3(1,0,0), Point3(0,1,0), T, inv(T))
Triangle(p1::Point3, p2::Point3, p3::Point3, T::Transformation) = Triangle(1, p1, p2, p3, T, inv(T))
Triangle(id::Integer, p1::Point3, p2::Point3, p3::Point3, T::Transformation) = Triangle(id, p1, p2, p3, T, inv(T))
Triangle(p1::Point3, p2::Point3, p3::Point3) = Triangle(p1, p2, p3, Transformation())


can_intersect(T::Triangle) = true
area(T::Triangle) = 0.5 * norm(cross(T.p2 - T.p1, T.p3 - T.p1))

function obj_bounds(T::Triangle)
    return BoundingBox(T.p1, T.p2, T.p3)
end

function world_bounds(T::Triangle)
    T.obj_to_world(obj_bounds(T))
end

function (T::Transformation)(C::Triangle)
    T2 = T * C.obj_to_world
    Triangle(C.id, C.p1, C.p2, C.p3, T2, inv(T2))
end


function shape_intersect(R::Ray, T::Triangle)
    ray = T.world_to_obj(R)
    
    e1 = T.p2 - T.p1
    e2 = T.p3 - T.p1
    
    s1 = cross(ray.direction, e2)
    divisor = dot(s1, e1)
    if divisor == 0.0
        return nothing, NaN, NaN
    end
    invd = 1.0 / divisor
    
    d = ray.origin - T.p1
    b1 = dot(d, s1) * invd
    if b1 < 0.0 || b1 > 1.0
        return nothing, NaN, NaN
    end
    
    s2 = cross(d, e1)
    b2 = dot(ray.direction, s2) * invd
    if b2 < 0.0 || b1+b2 > 1.0
        return nothing, NaN, NaN
    end
    
    t = dot(e2, s2) * invd
    if t < ray.tmin || t > ray.tmax
        return nothing, NaN, NaN
    end
    
    P = ray(t)
    
    n = normalize(Normal3(cross(e1, e2)))
    s = normalize(e1)
    
    n = flip_normal(n, ray.direction)
    
    DG = DifferentialGeometry(T.obj_to_world(P), T.obj_to_world(n), T.obj_to_world(s))
    
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




