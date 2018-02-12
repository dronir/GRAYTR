

# Interface of a Shape:
#  can_intersect(s::Shape), returns Bool 
#  area(s::Shape), returns Float64
#  obj_bounds(s::Shape), returns BoundingBox
#  world_bounds(s::Shape), returns BoundingBox
#  intersect(r::Ray, s::Shape) returns (Nullable{DifferentialGeometry}, tMin, rayEps)
#  intersectP(r::Ray, s::Shape) returns Bool

import Base.intersect

# Solve quadratic equation. Needed by different shapes' intersection calculation.
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

# Compute partial derivatives of normal vectors from Weingarten equations,
# given first and second derivatives of position vector.
# Used by several shape models.
function normal_derivatives(dpdu, dpdv, dpduu, dpduv, dpdvv)
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
    return dndu, dndv
end

get_shading_geometry(S::Shape, dg::DifferentialGeometry, T::Transformation) = dg

include("shapes/sphere.jl")
include("shapes/cylinder.jl")
include("shapes/disk.jl")
include("shapes/plane.jl")
include("shapes/box.jl")



