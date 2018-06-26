
# Sphere type and related methods

struct Sphere <: Shape
    id::Int64
    radius::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Sphere() = Sphere(1, 1.0, false, Transformation(), Transformation())
Sphere(r::Real, T::Transformation) = Sphere(1, r, false, T, inv(T))
Sphere(id::Int64, r::Real, T::Transformation) = Sphere(id, r, false, T, inv(T))
Sphere(id::Int64, r::Real, inverted::Bool, T::Transformation) = Sphere(id, r, inverted, T, inv(T))

can_intersect(s::Sphere) = true

area(s::Sphere) = 4π * s.radius^2

function obj_bounds(s::Sphere)
    v = Point3(s.radius, s.radius, s.radius)
    return BoundingBox(-v, v)
end

function world_bounds(s::Sphere)
    s.obj_to_world(obj_bounds(s))
end




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
    
    # Local coordinates
    phi = atan2(P.y, P.x)
    u = phi >= 0 ? phi / 2π : phi/2π + 1.0
    v = acos(clamp(P.z / sph.radius, -1.0, 1.0)) / π
    
    # Position vector partial derivatives
    zradius = sqrt(P.x^2 + P.y^2)
    cosphi = zradius ≈ 0.0 ? 1.0 : P.x / zradius
    sinphi = zradius ≈ 0.0 ? 0.0 : P.y / zradius
    dpdu = -Vector3(-2π*P.y, 2π*P.x, 0.0)
    dpdv = π * Vector3(P.z * cosphi, P.z * sinphi, -sph.radius * sin(v*π))
    
    # Normal vector partial derivatives through second derivatives of point
    dpduu = -4π^2 * Vector3(P.x, P.y, 0.0)
    dpduv = 2pi^2 * P.z * Vector3(-sinphi, cosphi, 0.0)
    dpdvv = -π^2 * P
    dndu, dndv = normal_derivatives(dpdu, dpdv, dpduu, dpduv, dpdvv)
    
    DG = DifferentialGeometry(
        sph.obj_to_world(P),
        u, v, sph,
        sph.obj_to_world(dpdu),
        sph.obj_to_world(dpdv),
        sph.obj_to_world(dndu),
        sph.obj_to_world(dndv)
    )
    return DG, t, 5e-4 * t
end

# Quick intersect function
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
