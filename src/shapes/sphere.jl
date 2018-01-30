
# Sphere type and related methods

struct Sphere <: Shape
    id::Int64
    radius::Float64
    inverted::Bool
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Sphere(id::Int64, r::Real, T::Transformation) = Sphere(id, r, false, T, inv(T))
Sphere(id::Int64, r::Real, inverted::Bool, T::Transformation) = Sphere(id, r, inverted, T, inv(T))

can_intersect(s::Sphere) = true

area(s::Sphere) = 4π * s.radius^2

function object_bounds(s::Sphere)
    v = Point3(s.radius, s.radius, s.radius)
    return BoundingBox(-v, v)
end

function world_bounds(s::Sphere)
    s.obj_to_world(object_bounds(s))
end

# Solve quadratic equation
function quadratic(A, B, C)
    dd = B^2 - 4*A*C
    if dd < 0.0
        return false, Inf, Inf
    end
    q = B<0 ? -0.5 * (B - sqrt(dd)) : -0.5 * (B + sqrt(dd))
    t0 = q/A
    t1 = C/q
    return true, min(t0,t1), max(t0,t1)
end


function intersect(r::Ray, sph::Sphere)
    ray = sph.world_to_obj(r)
    A = dot(ray.direction, ray.direction)
    B = 2.0 * dot(ray.direction, ray.origin)
    C = dot(ray.origin, ray.origin) - sph.radius^2
    hit, tnear, tfar = quadratic(A, B, C)
    
    # Not hit or both hits out of bounds
    if !hit || tnear > r.tmax || tfar < r.tmin || (tnear < r.tmin && tfar > r.tmax)
        return Nullable{DifferentialGeometry}(), NaN, NaN
    end
    t = tnear > r.tmin ? tnear : tfar
    P = ray(t)
    
    # Local coordinates
    X = P.x == 0.0 && P.y == 0.0 ? 1e-10*sph.radius : P.x
    phi = atan2(P.y, X)
    u = phi >= 0 ? phi / 2π : phi/2π + 1.0
    v = acos(clamp(P.z / sph.radius, -1.0, 1.0)) / π
    
    # Position vector partial derivatives
    zradius = sqrt(X^2 + P.y^2)
    cosphi = X / zradius
    sinphi = P.y / zradius
    dpdu = Vector3(-2π*P.y, 2π*X, 0.0)
    dpdv = π * Vector3(P.z * cosphi, P.z * sinphi, -sph.radius * sin(v*π))
    
    # Normal vector partial derivatives through second derivatives of point
    dpduu = -4π^2 * Vector3(X, P.y, 0.0)
    dpduv = 2pi^2 * P.z * Vector3(-sinphi, cosphi, 0.0)
    dpdvv = -π^2 * P
    E = dot(dpdu, dpdu)
    F = dot(dpdu, dpdv)
    G = dot(dpdv, dpdv)
    N = normalize(cross(dpdu, dpdv))
    e = dot(N, dpduu)
    f = dot(N, dpduv)
    g = dot(N, dpdvv)
    
    invEG = 1.0 / (E*G - F^2)
    dndu = Normal3((f*F - e*G) * invEG * dpdu + (e*F - f*E) * invEG * dpdv)
    dndv = Normal3((g*F - f*G) * invEG * dpdu + (f*F - g*E) * invEG * dpdv)
    
    
    DG = Nullable(DifferentialGeometry(
        sph.obj_to_world(P),
        sph.obj_to_world(N),
        u, v, sph,
        sph.obj_to_world(dpdu),
        sph.obj_to_world(dpdu),
        sph.obj_to_world(dndu),
        sph.obj_to_world(dndu)
    ))
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
