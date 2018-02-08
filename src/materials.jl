
include("brdf.jl")


struct BSDF
    dgs::DifferentialGeometry
    ng::Normal3
    nn::Normal3
    sn::Vector3
    tn::Vector3
    BxDFlist::Array{BxDF,1}
    refidx::Float64
end


function BSDF(dgs::DifferentialGeometry, ng::Normal3, L::Array{<:BxDF,1}, e::Float64)
    nn = dgs.n
    sn = normalize(dgs.dpdu)
    tn = cross(nn,sn)
    return BSDF(dgs, ng, nn, sn, tn, L, e)
end

BSDF(dgs::DifferentialGeometry, ng::Normal3, L::Array{BxDF,1}) = BSDF(dgs,ng,L,1.0)

ncomponents(B::BSDF) = size(B.BxDFlist)
world_to_local(B::BSDF, v::Vector3) = Vector3(dot(v, B.sn), dot(v, B.tn), dot(v, B.nn))
local_to_world(B::BSDF, v::Vector3) = Vector3(sn.x*v.x + tn.x*v.y + nn.x*vz,
                                              sn.y*v.x + tn.y*v.y + nn.y*vz,
                                              sn.z*v.x + tn.z*v.y + nn.z*vz)

function evaluate(B::BSDF, w0::Vector3, w1::Vector3)
    w0l = world_to_local(B, w0)
    w1l = world_to_local(B, w1)
    flags = 0
    if dot(w0, B.ng) * dot(w1, B.ng) > 0.0
        flags = flags & ~BSDF_TRANSMISSION
    else
        flags = flags & ~BSDF_REFLECTION
    end
    return sum(evaluate(some_BRDF, w0l, w1l) for some_BRDF in B.BxDFlist)
end


# The simplest material, just Lambertian

struct MatteMaterial <: Material
    Kd::Spectrum
end

function get_BSDF(M::MatteMaterial, dg::DifferentialGeometry, dgs::DifferentialGeometry)
    return BSDF(dgs, dg.n, [Lambert(M.Kd)], 1.0)
end
