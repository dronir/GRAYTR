
export AshikhminShirleySingle, Lambert, LommelSeeliger

# Utility functions
costheta(w::Vector3) = w.z
sintheta2(w::Vector3) = max(0.0, 1.0 - costheta(w)^2)
sintheta(w::Vector3) = sqrt(sintheta2(w))
cosphi(w::Vector3) = costheta(w) ≈ 1.0 ? 1.0 : clamp(w.x / sintheta(w), -1.0, 1.0)
sinphi(w::Vector3) = costheta(w) ≈ 1.0 ? 0.0 : clamp(w.y / sintheta(w), -1.0, 1.0)


"""
    evaluate(B::BxDF, world_to_local::Transformation, w0::Vector3, w1::Vector3)

Evaluate the value of a BRDF given a local transformation and the incident and emergent
directions. The local transformation is one that rotates the surface normal to coincide
with the +Z direction.

"""
function evaluate(B::BxDF, world_to_local::Transformation, w0::Vector3, w1::Vector3) 
    return evaluate(B, world_to_local(w0), world_to_local(w1))
end
    


######################################################################

"""
    Lambert{T<:Spectrum} <: BxDF

The Lambertian BRDF with spectrum `R`.

"""
struct Lambert{T<:Spectrum} <: BxDF
    R::T
end


"""
    evaluate(B::Lambert, w0::Vector3, w1::Vector3)

Evaluate the Lambertian BRDF for given direction vectors.

"""
evaluate(B::Lambert, w0::Vector3, w1::Vector3) = B.R ./ π


""""""
rho(B::Lambert, w0::Vector3) = B.R


"""
    generate_ray(B::Lambert, w0::Vector3) = nothing

Generate a new ray, weighed by the emergent light distribution of the Lambertian BRDF.

TODO: not implemented yet

"""
generate_ray(B::Lambert, w0::Vector3) = nothing


"""
    compute_pressure(B::Lambert, w0::Vector3, S::Spectrum)

Compute the radiation pressure on a surface

"""
function compute_pressure(B::Lambert, w0::Vector3, S::Spectrum)
    if costheta(w0) < 0.0
        return Vector3(0)
    end
    total_energy = integrate(B.R .* S)
    return -total_energy * (w0 + 2/3 * Z_AXIS)
end






######################################################################



"""
    LommelSeeliger{T<:Spectrum} <: BxDF

Lommel-Seeliger BRDF. Contains the surface spectrum and a phase function. The phase
function needs to be a one-parameter function defined on [0,π] and it's assumed to be
normalized so that its integral over the sphere equals 4π.

"""
struct LommelSeeliger{T<:Spectrum} <: BxDF
    R::T
    P::Function
end



"""
    evaluate(B::LommelSeeliger, w0::Vector3, w1::Vector3) 

Evaluate Lommel-Seeliger BRDF for given direction vectors.

TODO: check physics

"""
function evaluate(B::LommelSeeliger, w0::Vector3, w1::Vector3) 
    cos_alpha = clamp(dot(w0, w1), -1.0, 1.0)
    alpha = acos(cos_alpha)
    phase = B.P(alpha)
    return 0.25 * phase * B.R * max(0.0, 1.0 / (costheta(w0) + costheta(w1)))
end


SpectrumLike = Union{T,R} where T <: Spectrum where R <: Real


######################################################################


"""
    AshikhminShirleySingle{T<:Spectrum} <: BxDF

The Ashikhmin-Shirley BRDF where `R` is a spectrum, `d` is the diffuse component weigth
and `n` is the specular reflection width parameter. See Wetterer (2014).

"""
struct AshikhminShirleySingle{S<:SpectrumLike, T<:SpectrumLike} <: BxDF
    R::T
    specularity::S
    peak::Float64
end



"""
    Fresnel(R::Spectrum, costheta::Real, s::Real)

The Fresnel term of the Ashikhmin-Shirley BRDF (see Wetterer, 2014).
"""
@inline Fresnel(Rs::SpectrumLike, x::Real) = Rs .+ (1.0 .- Rs) .* (1.0 - x)^5


"""
    BlinnPhong(n::Real, hdn::Real)
    
The Blinn-Phong function used in the Ashikhmin-Shirley BRDF (see Wetterer, 2014).
"""
@inline BlinnPhong(n::SpectrumLike, hdn::Real) = (n .+ 1) .* hdn^n ./ 8π



@inline f5(x) = (one(x) - (one(x) - x/2)^5)


"""
    evaluate(B::AshikhminShirleySingle, w0::Vector3, w1::Vector3)

Evaluate Ashikhmin-Shirley BRDF for given direction vectors.

"""
function evaluate(B::AshikhminShirleySingle, w0::Vector3, w1::Vector3)
    h = normalize((w0 + w1) / 2.0)
    D = BlinnPhong(B.peak, h.z)
    F = Fresnel(B.R, costheta(w1))
    specular = F .* (D ./ (dot(h,w1) * max(costheta(w0), costheta(w1))))
    diffuse = B.R .* (1 .- B.specularity) .* (28/23π * f5(costheta(w0)) * f5(costheta(w1)))
    
    return diffuse .+ specular
end




######################################################################


"""
    SpecularDiffuse{T<:Spectrum} <: BxDF

Specular + Diffuse reflection. The fraction of diffusely reflected light is `fd`, the
fraction of specularly reflected light is `fs`. Absorption is `1 - fd - fs`.

TODO: Make all the parameters spectral.

"""
struct SpecularDiffuse{T<:Spectrum} <: BxDF
    R::T
    fd::Float64
    fs::Float64
end



"""
    evaluate(B::SpecularDiffuse, w0::Vector3, w1::Vector3)

Evaluate diffuse + specular BRDF for given direction vectors.

TODO: normalization is probably off...

"""
function evaluate(B::SpecularDiffuse, w0::Vector3, w1::Vector3)
    if normalize(w0 + w1) ≈ Z_AXIS
        return (B.fd / π + B.fs) .* B.R
    else
        return B.fd / π .* B.R
    end
end


"""
    compute_pressure(B::SpecularDiffuse, w0::Vector3)

Compute the radiation pressure on a surface.

"""
function compute_pressure(B::SpecularDiffuse, w0::Vector3, S::Spectrum)
    if costheta(w0) < 0.0
        return Vector3(0)
    end
    fa = 1.0 - B.fd - B.fs
    total_energy = integrate(B.R .* S)
    return -total_energy * ((fa + B.fd) * (w0 + 2/3 * Z_AXIS) + 2*B.fs*costheta(w0) * Z_AXIS)
end
