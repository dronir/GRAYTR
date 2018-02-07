import Base.+, Base.*, Base./

struct NoLight <: Spectrum end

+(N::NoLight, S::Spectrum) = S
+(S::Spectrum, N::NoLight) = S
*(N::NoLight, S::Spectrum) = S
*(S::Spectrum, N::NoLight) = S
/(N::NoLight, S::Spectrum) = S
/(S::Spectrum, N::NoLight) = S


struct SampledSpectrum <: Spectrum
    low::Float64
    high::Float64
    values::Array{Float64,1}
end


function +(S1::SampledSpectrum, S2::SampledSpectrum)
    if !(S1.low ≈ S2.low && S1.high ≈ S2.high) 
        error("Spectrum sizes don't match.")
    end
    return SampledSpectrum(S1.values + S2.values)
end
+(S::SampledSpectrum, c::Real) = SampledSpectrum(S.low, S.high, S.values + c)
+(c::Real, S::SampledSpectrum) = S+c

isblack(S::SampledSpectrum) = max(SampledSpectrum.values) ≈ 0.0

import Base.broadcast
broadcast(f, S::SampledSpectrum) = SampledSpectrum(S.low, S.high, f.(S.values))


function XYZtoRGB(xyz::Vector3)
    return Vector3(
        3.240479*xyz.x - 1.537150*xyz.y - 0.498535*xyz.z,
       -0.969256*xyz.x + 1.875991*xyz.y + 0.041556*xyz.z,
        0.055648*xyz.x - 0.204043*xyz.y + 1.057311*xyz.z)
end

function RGBtoXYZ(rgb::Vector3)
    return Vector3(
        0.412453*rgb.x + 0.357580*rgb.y + 0.180423*rgb.z,
        0.212671*rgb.x + 0.715160*rgb.y + 0.072169*rgb.z,
        0.019334*rgb.x + 0.119193*rgb.y + 0.950227*rgb.z)
end


struct RGBSpectrum <: Spectrum
    r::Float64
    g::Float64
    b::Float64
end

isblack(S::RGBSpectrum) = S.r == 0.0 && S.g == 0.0 && S.b == 0.0

to_XYZ(S::RGBSpectrum) = [x for x in (Vector3(S.r, S.g, S.b))]

+(a::RGBSpectrum, b::RGBSpectrum) = RGBSpectrum(a.r+b.r, a.g+b.g, a.b+b.b)
*(x::Number, S::RGBSpectrum) = RGBSpectrum(x*S.r, x*S.g, x*S.b)
*(S::RGBSpectrum, x::Number) = x*S
/(S::RGBSpectrum, x::Number) = RGBSpectrum(S.r/x, S.g/x, S.b/x)
*(a::RGBSpectrum, b::RGBSpectrum) = RGBSpectrum(a.r*b.r, a.g*b.g, a.b*b.b)
broadcast(f, S::RGBSpectrum) = RGBSpectrum(f(S.r), f(S.g), f(S.b))
