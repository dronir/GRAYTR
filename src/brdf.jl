
# Utility functions
costheta(w::Vector3) = w.z
sintheta2(w::Vector3) = max(0.0, 1.0 - costheta(w)^2)
sintheta(w::Vector3) = sqrt(sintheta2(w))
cosphi(w::Vector3) = sintheta(w) ≈ 0.0 ? 1.0 : clamp(w.x / sintheta(w), -1.0, 1.0)
sinphi(w::Vector3) = sintheta(w) ≈ 0.0 ? 0.0 : clamp(w.y / sintheta(w), -1.0, 1.0)




# Interface for BxDF subtypes:
#  BSDF_type(B::BxDF) returns a type (TODO: how to parametrize?)
#  f(B::BxDF, w0::Vector3, w1::Vector3) returns the value
#  sample_f(B::BxDF)
#  rho(B::BxDF, w0::Vector3) returns Spectrum (hemispherical reflectance)

const BSDF_REFLECTION = 1
const BSDF_TRANSMISSION = 2
const BSDF_SPECULAR = 4
const BSDF_DIFFUSE = 8
const BSDF_GLOSSY = 16
const BSDF_ALL_TYPES = BSDF_DIFFUSE | BSDF_GLOSSY | BSDF_SPECULAR
const BSDF_ALL_REFLECTION = BSDF_REFLECTION | BSDF_ALL_TYPES
const BSDF_ALL_TRANSMISSION = BSDF_TRANSMISSION | BSDF_ALL_TYPES
const BSDF_ALL = BSDF_ALL_TRANSMISSION | BSDF_ALL_REFLECTION


# Lambertian

struct Lambert{T<:Spectrum} <: BxDF
    R::T
end

BSDF_type(B::Lambert) = BSDF_REFLECTION | BSDF_DIFFUSE
evaluate(B::Lambert, w0::Vector3, w1::Vector3) = B.R / π
rho(B::Lambert, w0::Vector3) = B.R



# Lommel-Seeliger

struct LommelSeeliger{T<:Spectrum} <: BxDF
    R::T
end

BSDF_type(B::LommelSeeliger) = BSDF_REFLECTION | BSDF_DIFFUSE
evaluate(B::LommelSeeliger, w0::Vector3, w1::Vector3) = B.R .* max(0.0, 1.0 / (costheta(w0) + costheta(w1)))



# Ashkhmin-Shirley

struct AshkhminShirleySingle{T<:Spectrum} <: BxDF
    R::T
    n::Float64
end


BSDF_type(B::AshkhminShirleySingle) = BSDF_REFLECTION | BSDF_DIFFUSE

#Fresnel(R::Spectrum, costheta::Real) = R + (1 - R) * (1 - costheta)^5
Fresnel(R::Spectrum, costheta::Real) = R*(1 - (1 - costheta)^5) + (1 - costheta)^5
BlinnPhong(e::Real, hdn::Real) = (e+2) * hdn^e / 2π

function evaluate(B::AshkhminShirleySingle, w0::Vector3, w1::Vector3)
    h = normalize((w0 + w1) / 2.0)
    D = BlinnPhong(B.n, h.z)
    specular = Fresnel(B.R, costheta(w1)) * (D / (4π * dot(h,w1) * max(costheta(w0), costheta(w1))))
    diffuse = (28/23π * (1 - (1 - costheta(w0)/2)^5) * (1 - (1 - costheta(w1)/2)^5)) .* (B.R - B.R.^2)
    return specular + diffuse
end

