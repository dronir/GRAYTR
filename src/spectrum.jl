import Base.+, Base.*, Base./


# ------------------------------------------------
# Load CIE spectra

raw_CIE_data = readcsv("lin2012xyz2e_1_7sf.csv")
const CIE_LAMBDA = raw_CIE_data[:,1]
const CIE_X = raw_CIE_data[:,2]
const CIE_Y = raw_CIE_data[:,3]
const CIE_Z = raw_CIE_data[:,4]


# ------------------------------------------------
# The NoLight dummy spectrum

struct NoLight <: Spectrum end

isblack(N::NoLight) = true

*(x::Number, N::NoLight) = N
+(N::NoLight, M::NoLight) = NoLight()
+(N::NoLight, S::Spectrum) = S
+(S::Spectrum, N::NoLight) = S
*(N::NoLight, S::Spectrum) = NoLight()
*(S::Spectrum, N::NoLight) = NoLight()
/(N::NoLight, S::Spectrum) = S
/(S::Spectrum, N::NoLight) = S


# ------------------------------------------------
# The SampledSpectrum type

struct SampledSpectrum <: Spectrum
    low::Int64
    high::Int64
    values::Array{Float64,1}
end

isblack(S::SampledSpectrum) = maximum(S.values) ≈ 0.0

# Addition of SampledSpectrum objects
function +(S1::SampledSpectrum, S2::SampledSpectrum)
    if !(S1.low == S2.low && S1.high == S2.high && size(S1.values) == size(S2.values)) 
        error("Spectrum sizes don't match.")
    end
    return SampledSpectrum(S1.values + S2.values)
end
+(S::SampledSpectrum, c::Real) = SampledSpectrum(S.low, S.high, S.values + c)
+(c::Real, S::SampledSpectrum) = S+c


# Multiplication of SampledSpectrum objects
function *(S1::SampledSpectrum, S2::SampledSpectrum)
    if !(S1.low == S2.low && S1.high == S2.high && size(S1.values) == size(S2.values))
        error("Spectrum sizes don't match.")
    end
    return SampledSpectrum(S1.low, S1.high, S1.values .* S2.values)
end
*(S::SampledSpectrum, c::Real) = c*S
*(c::Real, S::SampledSpectrum) = SampledSpectrum(S.low, S.high, c*S.values)


# Various functions
import Base.broadcast, Base.size, Base.length
broadcast(f, S::SampledSpectrum) = SampledSpectrum(S.low, S.high, f.(S.values))
size(S::SampledSpectrum) = size(S.values)
length(S::SampledSpectrum) = length(S.values)

function interpolate(S::SampledSpectrum, a::Real)
    if a < S.low
        return S.values[1]
    elseif a >= S.high
        return S.values[end]
    end
    delta = (S.high - S.low) / (length(S.values) - 1)
    N = Int64(fld(a - S.low, delta))
    return lerp((a - S.low - N*delta) / delta, S.values[N+1], S.values[N+2])
end

function interpolate(S::SampledSpectrum, a::Real, b::Real)
    if b < S.low
        return S.values[1]
    elseif a > S.high
        return S.values[end]
    end
    
    
end


# TODO: Convert SampledSpectrum to XYZ
#function to_XYZ(S::SampledSpectrum)
#end


# ------------------------------------------------
# Conversion between XYZ and RGB

function XYZtoRGB(xyz::Array{Float64,1})
    return [3.240479*xyz[1] - 1.537150*xyz[2] - 0.498535*xyz[3],
           -0.969256*xyz[1] + 1.875991*xyz[2] + 0.041556*xyz[3],
            0.055648*xyz[1] - 0.204043*xyz[2] + 1.057311*xyz[3]]
end  

function RGBtoXYZ(rgb::Array{Float64,1})
    return [0.412453*rgb[1] + 0.357580*rgb[2] + 0.180423*rgb[3],
            0.212671*rgb[1] + 0.715160*rgb[2] + 0.072169*rgb[3],
            0.019334*rgb[1] + 0.119193*rgb[2] + 0.950227*rgb[3]]
end


# ------------------------------------------------
# The RGBSpectrum type and 

struct RGBSpectrum <: Spectrum
    r::Float64
    g::Float64
    b::Float64
end

isblack(S::RGBSpectrum) = S.r == 0.0 && S.g == 0.0 && S.b == 0.0

to_XYZ(S::RGBSpectrum) = RGBtoXYZ([S.r, S.g, S.b])

+(a::RGBSpectrum, b::RGBSpectrum) = RGBSpectrum(a.r+b.r, a.g+b.g, a.b+b.b)
*(x::Number, S::RGBSpectrum) = RGBSpectrum(x*S.r, x*S.g, x*S.b)
*(S::RGBSpectrum, x::Number) = x*S
/(S::RGBSpectrum, x::Number) = RGBSpectrum(S.r/x, S.g/x, S.b/x)
*(a::RGBSpectrum, b::RGBSpectrum) = RGBSpectrum(a.r*b.r, a.g*b.g, a.b*b.b)
broadcast(f, S::RGBSpectrum) = RGBSpectrum(f(S.r), f(S.g), f(S.b))
