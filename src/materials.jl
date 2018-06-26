
include("brdf.jl")


function local_transformation(dg::DifferentialGeometry)
    # Make an orthonormal basis of the normal vector, a vector tangent to the surface,
    # and their cross product.
    n = dg.n
    s = normalize(dg.dpdu)
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



function evaluate(B::BxDF, world_to_local::Transformation, w0::Vector3, w1::Vector3)
    w0l = world_to_local(w0)
    w1l = world_to_local(w1)
    flags = 0
    if sign(w0l.z) * sign(w1l.z) < 0.0
        flags = flags & ~BSDF_TRANSMISSION
    else
        flags = flags & ~BSDF_REFLECTION
    end
    return evaluate(B, w0l, w1l)
end

