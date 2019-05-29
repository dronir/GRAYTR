
struct Disk <: Shape
    id::Int64
    rmin::Float64
    rmax::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Disk(T::Transformation) = Disk(1, 0.0, 1.0, false, T, inv(T))
Disk() = Disk(1, 0.0, 1.0, false, Transformation(), Transformation())

can_intersect(D::Disk) = true
area(D::Disk) = π * (D.rmax^2 - D.rmin^2)

function obj_bounds(D::Disk)
    pMin = Point3(-D.rmax, -D.rmax, 0.0)
    pMax = Point3( D.rmax,  D.rmax, 0.0)
    return BoundingBox(pMin, pMax)
end

function world_bounds(D::Disk)
    D.obj_to_world(obj_bounds(D))
end

function (T::Transformation)(D::Disk)
    T2 = T * D.obj_to_world
    Disk(D.id, D.rmin, D.rmax, D.inverted, T2, inv(T2))
end



function shape_intersect(R::Ray, D::Disk)
    ray = D.world_to_obj(R)
    if ray.direction.z ≈ 0.0
        return nothing, NaN, NaN
    end
    t = -ray.origin.z / ray.direction.z
    if t < ray.tmin || t > ray.tmax
        return nothing, NaN, NaN 
    end
    P = ray(t)
    r2 = P.x^2 + P.y^2
    if r2 > D.rmax^2 || r2 < D.rmin^2
        return nothing, NaN, NaN
    end

    dpdu = Vector3(-P.y, P.x, 0.0)

    n =Normal3(0, 0, sign(ray.origin.z))
    s =normalize(dpdu)

    DG = DifferentialGeometry(D.obj_to_world(P), D.obj_to_world(n), D.obj_to_world(s))
    
    return DG, t, 5e-4 * t
end

function intersectP(R::Ray, D::Disk)
    ray = D.world_to_obj(R)
    if ray.direction.z ≈ 0.0
        return false
    end
    t = -ray.origin.z / ray.direction.z
    if t < ray.tmin || t > ray.tmax
        return false
    end
    P = ray(t)
    r2 = P.x^2 + P.y^2
    if r2 > D.rmax^2 || r2 < D.rmin^2
        return false
    end
    return true
end

