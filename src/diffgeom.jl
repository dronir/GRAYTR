
# Differential geometry.

# This struct represents the geometry on a point on a surface where a ray has hit.

struct DifferentialGeometry
    # Location and normal vector on the surface
    p::Point3
    n::Normal3
    T::Transformation
end


function DifferentialGeometry(p::Point3, n::Normal3, s::Vector3)
    n = normalize(n)
    s = normalize(s)
    return DifferentialGeometry(p, n, local_transformation(n, s))
end




"""
    flip_normal(n::Normal3, d::Vector3)

Inverts the normal vector `n` if the normal and the direction vector `d` are in the same
hemisphere. This is to always get normals pointing "out" of a surface regardless of whether
the ray comes from "inside" or "outside" the shape.

"""
flip_normal(n::Normal3, d::Vector3) = dot(n,d) < 0.0 ? n : -n



"""
    local_transformation(n::Normal3, s::Vector3)

Get the transformation which transforms a vector from the world coordinates to
surface-local coordinates, where the Z-axis is along the surface normal.

Makes an orthonormal basis of the normal vector `n`, a vector `s` tangent to the surface,
and their cross product. Both vectors are assumed to be normalized.

"""
function local_transformation(n::Normal3, s::Vector3)
    t = cross(n, s)
    M = zeros(4,4)
    for i = 1:3
        M[1,i] = s[i]
        M[2,i] = t[i]
        M[3,i] = n[i]
    end
    M[4,4] = 1.0
    
    return Transformation(M, M')
end

