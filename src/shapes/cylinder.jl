
struct Cylinder <: Shape
    id::Int64
    radius::Float64
    zmin::Float64
    zmax::Float64
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Cylinder(T::Transformation) = Cylinder(1, 1.0, -1.0, 1.0, T, inv(T))
Cylinder(id::Int64, radius::Real, zmin::Real, zmax::Real) = Cylinder(id, radius, zmin, zmax, Transformation(), Transformation())


can_intersect(C::Cylinder) = true
area(C::Cylinder) = 2π * C.radius * (C.zmax - C.zmin)

function obj_bounds(C::Cylinder)
    pMin = Point3(-C.radius, -C.radius, C.zmin)
    pMax = Point3( C.radius,  C.radius, C.zmax)
    return BoundingBox(pMin, pMax)
end

function world_bounds(C::Cylinder)
    C.obj_to_world(obj_bounds(C))
end

function (T::Transformation)(C::Cylinder)
    T2 = T * C.obj_to_world
    Cylinder(C.id, C.radius, C.zmin, C.zmax, T2, inv(T2))
end


function cylinder_coefs(ray::Ray, cyl::Cylinder)
    A = ray.direction.x^2 + ray.direction.y^2
    B = 2.0 * (ray.direction.x * ray.origin.x + ray.direction.y * ray.origin.y)
    C = ray.origin.x^2 + ray.origin.y^2 - cyl.radius^2
    return A, B, C
end


function shape_intersect(R::Ray, cyl::Cylinder)
    ray = cyl.world_to_obj(R)
    
    A, B, C = cylinder_coefs(ray, cyl)
    if A ≈ 0.0
        return nothing, NaN, NaN
    end
    hit, tnear, tfar = quadratic(A,B,C)
    
    if !hit || tnear > ray.tmax || tfar < ray.tmin || (tnear < ray.tmin && tfar > ray.tmax)
        # The ray didn't even hit the infinitely long extension of the cylinder.
        return nothing, NaN, NaN
    end
    
    P1 = ray(tnear)
    P2 = ray(tfar)
    
    if cyl.zmin <= P1.z <= cyl.zmax && tnear >= ray.tmin
        P = P1
        t = tnear
    elseif cyl.zmin <= P2.z <= cyl.zmax && tnear <= ray.tmax
        P = P2
        t = tfar
    else
        return nothing, NaN, NaN
    end
        
    n = normalize(Normal3(P.x, P.y, 0.0))
    s = Vector3(0.0, 0.0, 1.0)
    
    n = flip_normal(n, ray.direction)
    
    DG = DifferentialGeometry(cyl.obj_to_world(P), cyl.obj_to_world(n), cyl.obj_to_world(s))
    
    return DG, t, 5e-4 * t
end

function intersectP(r::Ray, cyl::Cylinder)
    ray = cyl.world_to_obj(r)
    A, B, C = cylinder_coefs(ray, cyl)
    if A ≈ 0.0
        return false
    end
    hit, tnear, tfar = quadratic(A,B,C)
    
    if !hit || tnear > ray.tmax || tfar < ray.tmin || (tnear < ray.tmin && tfar > ray.tmax)
        return false
    end
    
    P1 = ray(tnear)
    P2 = ray(tfar)
    
    return (cyl.zmin <= P1.z <= cyl.zmax && tnear >= ray.tmin) || (cyl.zmin <= P2.z <= cyl.zmax && tnear <= ray.tmax)
        
        
end

