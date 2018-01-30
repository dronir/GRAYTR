
# Differential geometry.

# This struct represents the geometry on a point on a surface where a ray has hit.

struct DifferentialGeometry
    # Location and normal vector on the surface
    p::Point3
    n::Normal3
    # Surface coordinates
    u::Float64
    v::Float64
    # Reference to the Shape the point is on
    shape::Shape
    # Partial derivatives of the point and normal wrt surface coordinates
    dpdu::Vector3
    dpdv::Vector3
    dndu::Normal3
    dndv::Normal3
end

# Constructor that computes the normal vector.
function DifferentialGeometry(p::Point3, u::Float64, v::Float64,
                              dpdu::Vector3, dpdv::Vector3, 
                              dndu::Normal3, dndv::Normal3, sh::Shape)
    normal = normalize(cross(dndu, dndv))
    if sh.inverted ⊻ swaps_handedness(sh.obj_to_world)
        normal -= normal
    end
    DifferentialGeometry(p, normal, u, v, sh, dpdu, dpdv, dndu, dndv)
end

