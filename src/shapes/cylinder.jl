
struct Cylinder <: Shape
    id::Int64
    radius::Float64
    zmin::Float64
    zmax::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Cylinder(T::Transformation) = Cylinder(1, 1.0, 0.0, 1.0, false, T, inv(T))

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
    Cylinder(C.id, C.radius, C.zmin, C.zmax, C.inverted, T2, inv(T2))
end



function shape_intersect(R::Ray, cyl::Cylinder)
    ray = cyl.world_to_obj(R)
    A = ray.direction.x^2 + ray.direction.y^2
    if A ≈ 0.0
        return nothing, NaN, NaN
    end
    B = 2.0 * (ray.direction.x * ray.origin.x + ray.direction.y * ray.origin.y)
    C = ray.origin.x^2 + ray.origin.y^2 - cyl.radius^2
    hit, tnear, tfar = quadratic(A,B,C)
    
    if !hit || tnear > ray.tmax || tfar < ray.tmin || (tnear < ray.tmin && tfar > ray.tmax)
        # The ray didn't hit the infinitely long extension of the cylinder.
        # This it misses the cylinder itself.
        return nothing, NaN, NaN
    end
    
    t = tnear > ray.tmin ? tnear : tfar
    P = ray(t)
    
    if P.z < cyl.zmin || P.z > cyl.zmax
        # The z coordinate of the intersection with the infinite cylinder is outside
        # of the actual cylinder's extent.
        return nothing, NaN, NaN
    end
    
    phi = atan(P.y, P.x)
    phi = phi >= 0.0 ? phi : phi + 2π
    u = phi / 2π
    v = (P.z - cyl.zmin) / (cyl.zmax - cyl.zmin)
    
    dpdu = Vector3(-2π * P.y, 2π * P.x, 0.0)
    dpdv = Vector3(0.0, 0.0, cyl.zmax - cyl.zmin)
    dpduu = Vector3(-4π^2*P.x, -4π^2*P.y, 0.0)
    dpduv = Vector3(0.0)
    dpdvv = Vector3(0.0)
    dndu, dndv = normal_derivatives(dpdu, dpdv, dpduu, dpduv, dpdvv)
    
    DG = DifferentialGeometry(
        cyl.obj_to_world(P),
        u, v, cyl,
        cyl.obj_to_world(dpdu),
        cyl.obj_to_world(dpdv),
        cyl.obj_to_world(dndu),
        cyl.obj_to_world(dndv)
    )
    return DG, t, 5e-4 * t
end

function intersectP(r::Ray, cyl::Cylinder)
    ray = cyl.world_to_obj(r)
    A = ray.direction.x^2 + ray.direction.y^2
    B = 2.0 * (ray.direction.x * ray.origin.x + ray.direction.y * ray.origin.y)
    C = ray.origin.x^2 + ray.origin.y^2 - cyl.radius^2
    hit, tnear, tfar = quadratic(A,B,C)
    
    if !hit || tnear > ray.tmax || tfar < ray.tmin || (tnear < ray.tmin && tfar > ray.tmax)
        return false
    end
    
    t = tnear > ray.tmin ? tnear : tfar
    P = ray(t)
    
    return P.z >= cyl.zmin && P.z <= cyl.zmax
end

