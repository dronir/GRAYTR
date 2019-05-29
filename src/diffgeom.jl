
# Differential geometry.

# This struct represents the geometry on a point on a surface where a ray has hit.

struct DifferentialGeometry
    # Location and normal vector on the surface
    p::Point3
    n::Normal3
    # Surface coordinates
    u::Float64
    v::Float64

    # Partial derivatives of the point and normal wrt surface coordinates
    dpdu::Vector3
    dpdv::Vector3
    dndu::Normal3
    dndv::Normal3
end

# Constructor that computes the normal vector.
function DifferentialGeometry(p::Point3, u::Float64, v::Float64, sh::Shape,
                              dpdu::Vector3, dpdv::Vector3, 
                              dndu::Normal3, dndv::Normal3)
    normal = normalize(cross(dpdu, dpdv))
    if sh.inverted ⊻ swaps_handedness(sh.obj_to_world)
        normal = -normal
    end
    DifferentialGeometry(p, normal, u, v, dpdu, dpdv, dndu, dndv)
end




"""
    local_transformation(dg::DifferentialGeometry)

Get the transformation which transforms a vector from the world coordinates to
surface-local coordinates, where the Z-axis is along the surface normal.

"""
function local_transformation(dg::DifferentialGeometry)
    # Make an orthonormal basis of the normal vector, a vector tangent to the surface,
    # and their cross product.
    s = normalize(dg.dpdu)
    t = cross(dg.n, s)
    
    M = zeros(4,4)
    for i = 1:3
        M[1,i] = s[i]
        M[2,i] = t[i]
        M[3,i] = dg.n[i]
    end
    M[4,4] = 1.0
    
    return Transformation(M, M')
end

