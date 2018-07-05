import Base.+, Base.-, Base.*, Base./, Base.convert, Base.zero, Base.one


# ------------------------------------------------
# Load CIE spectra

raw_CIE_data = readcsv("lin2012xyz2e_1_7sf.csv")
const N_CIE = length(raw_CIE_data[:,1])
const CIE_LAMBDA = vec(raw_CIE_data[:,1])
const CIE_X = vec(raw_CIE_data[:,2])
const CIE_Y = vec(raw_CIE_data[:,3])
const CIE_Z = vec(raw_CIE_data[:,4])





# ------------------------------------------------
# The NoLight dummy spectrum

struct NoLight <: Spectrum end

const nolight = NoLight()

isblack(N::NoLight) = true
to_XYZ(S::NoLight) = [0.0, 0.0, 0.0]

*(x::Number, N::NoLight) = nolight
+(N::NoLight, M::NoLight) = nolight
+(N::NoLight, S::Spectrum) = S
+(S::Spectrum, N::NoLight) = S
#*(N::NoLight, S::Spectrum) = nolight
#*(S::Spectrum, N::NoLight) = nolight
#/(N::NoLight, S::Spectrum) = S


# ------------------------------------------------
# The SampledSpectrum type

struct SampledSpectrum{N} <: Spectrum
    low::Int64
    high::Int64
    values::Array{Float64,1}
    function SampledSpectrum{N}(low::Int64, high::Int64, values::Array{Float64,1}) where N 
        length(values) != N ? error("Type parameter N not equal to length(values): $N != $(length(values))") : new(low, high, values)
    end
end

one(S::SampledSpectrum{N}) where N = SampledSpectrum{N}(S.low, S.high, ones(Float64, N))
zero(S::SampledSpectrum{N}) where N = SampledSpectrum{N}(S.low, S.high, zeros(Float64, N))
convert(T::Type{SampledSpectrum{N}}, L::NoLight) where N = SampledSpectrum(200, 1000, zeros(Float64, N))

SampledSpectrum(low::Integer, high::Integer, values::Array{Float64,1}) = SampledSpectrum{length(values)}(low, high, values)

isblack(S::SampledSpectrum) = maximum(S.values) ≈ 0.0

# Addition of SampledSpectrum objects
function +(S1::SampledSpectrum{N}, S2::SampledSpectrum{N}) where N
    if !(S1.low == S2.low && S1.high == S2.high) 
        error("Spectrum limits don't match:\n$S1 \n $S2")
    end
    return typeof(S1)(S1.low, S2.high, S1.values + S2.values)
end
+(S::SampledSpectrum{N}, c::Real) where N = SampledSpectrum{N}(S.low, S.high, S.values + c)
+(c::Real, S::SampledSpectrum) = S+c
-(c::Real, S::SampledSpectrum{N}) where N = SampledSpectrum{N}(S.low, S.high, c - S.values)

function -(S1::SampledSpectrum{N}, S2::SampledSpectrum{N}) where N
    if !(S1.low == S2.low && S1.high == S2.high) 
        error("Spectrum limits don't match:\n$S1 \n $S2")
    end
    return typeof(S1)(S1.low, S2.high, S1.values - S2.values)
end


# Multiplication of SampledSpectrum objects with reals
function *(S1::SampledSpectrum{N}, S2::SampledSpectrum{N}) where N
    if !(S1.low == S2.low && S1.high == S2.high) 
        error("Spectrum limits don't match:\n$S1 \n $S2")
    end
    return typeof(S1)(S1.low, S2.high, S1.values .* S2.values)
end
*(S::SampledSpectrum, c::Real) = c*S
*(c::Real, S::SampledSpectrum{N}) where N = SampledSpectrum{N}(S.low, S.high, c*S.values)
/(S::SampledSpectrum{N}, c::Real) where N = SampledSpectrum{N}(S.low, S.high, S.values/c)

# Various functions
import Base.broadcast, Base.size, Base.length
#broadcast(f::Function, S::SampledSpectrum{N}) where N = SampledSpectrum{N}(S.low, S.high, broadcast(f, S.values))
#broadcast(f::Function, S::SampledSpectrum{N}, Y::SampledSpectrum{N}) where N = SampledSpectrum{N}(S.low, S.high, broadcast(f, S.values, Y.values))
#broadcast(f::Function, S::SampledSpectrum{N}, y...) where N = SampledSpectrum{N}(S.low, S.high, broadcast(f, S.values, y...))
size(S::SampledSpectrum{N}) where N = N
length(S::SampledSpectrum{N}) where N = N

function interpolate(S::SampledSpectrum, a::Real)
    if a < S.low
        return S.values[1]
    elseif a >= S.high
        return S.values[end]
    end
    delta = (S.high - S.low) / (length(S) - 1)
    N = Int64(fld(a - S.low, delta))
    return lerp((a - S.low - N*delta) / delta, S.values[N+1], S.values[N+2])
end

function interpolate(S::SampledSpectrum, a::Real, b::Real)
    # We are assuming that a < b
    if b <= S.low
        return S.values[1]
    elseif a >= S.high
        return S.values[end]
    end
    N = length(S)
    delta = (S.high - S.low) / (N - 1)

    # The N values tell how many full deltas are between S.low and a or b
    Na = Int64(fld(a - S.low, delta))
    Nb = Int64(fld(b - S.low, delta))

    # Cases:
    # 1. Both a and b under the low limit; handled above
    # 2. Both a and b over the high limit; handled above
    # 3. a under low, b inside the limits: "low" Ia + mid + Ib
    # 4. a inside the limits, b over high: Ia + mid + "high" Ib
    # 5. both inside the limits, same bin: special case
    # 6. both inside the limits, different bin:  Ia + mid + Ib
    # 7. a under low and b over high: "low" Ia + mid + "high" Ib
    
    # Case 5:
    if Na == Nb
        # a and b are in the same spectrum interval, just linearly integrate between them.
        ya = lerp((a - S.low - Na*delta) / delta, S.values[Na+1], S.values[Na+2])
        yb = lerp((b - S.low - Na*delta) / delta, S.values[Na+1], S.values[Na+2])
        return 0.5 * (ya + yb)
    end
    
    if Na < 0
        # Cases 3 and 7
        # a is below the low end, first integral is constant low value over [a, S.low]
        Ia = S.values[1] * (S.low - a)
        k0 = 1
    else
        # Cases 4 and 6
        da = (S.low + (Na+1)*delta) - a
        ya = lerp(1 - da / delta, S.values[Na+1], S.values[Na+2])
        Ia = 0.5 * da * (ya + S.values[Na+2])
        k0 = Na+2
    end
    if Nb >= length(S.values)-1
        # Cases 4 and 7
        # b is above high end, last integral is constant high value over [S.high, b]
        Ib = S.values[end] * (b - S.high)
        k1 = N-1
    else
        # Cases 3 and 6
        db = b - (S.low + Nb*delta)
        yb = lerp(db / delta, S.values[Nb+1], S.values[Nb+2])
        Ib = 0.5 * db * (S.values[Nb+1] + yb)
        k1 = Nb
    end
    
    Imid = 0.0
    if (Nb-Na) != 1
        for i = k0:k1
            Imid += 0.5*delta*(S.values[i] + S.values[i+1])
        end
    end
    
    return (Ia + Imid + Ib) / (b - a)
end


function to_XYZ(S::SampledSpectrum)
    out = [0.0, 0.0, 0.0]
    for i = 1:N_CIE
        lam = CIE_LAMBDA[i]
        s = interpolate(S, lam)
        out[1] += s * CIE_X[i]
        out[2] += s * CIE_Y[i]
        out[3] += s * CIE_Z[i]
    end
    return out
end


# ------------------------------------------------
# The RGBSpectrum type and 

struct RGBSpectrum <: Spectrum
    r::Float64
    g::Float64
    b::Float64
end

import Base.getindex
getindex(S::RGBSpectrum, i::Integer) = i==1 ? S.r : i==2 ? S.g : i==3 ? S.b : error("Out ouf bounds index for RGBSpectrum: $i")

zero(T::Type{RGBSpectrum}) = RGBSpectrum(0.0, 0.0, 0.0)
zero(T::RGBSpectrum) = RGBSpectrum(0.0, 0.0, 0.0)
convert(T::Type{RGBSpectrum}, N::NoLight) = RGBSpectrum(0.0, 0.0, 0.0)

isblack(S::RGBSpectrum) = S.r == 0.0 && S.g == 0.0 && S.b == 0.0

to_XYZ(S::RGBSpectrum) = RGBtoXYZ(S)
to_RGB(S::RGBSpectrum) = S

+(a::RGBSpectrum, b::RGBSpectrum) = RGBSpectrum(a.r+b.r, a.g+b.g, a.b+b.b)
-(n::Real, S::RGBSpectrum) = RGBSpectrum(1-S.r, 1-S.g, 1-S.b)
*(x::Number, S::RGBSpectrum) = RGBSpectrum(x*S.r, x*S.g, x*S.b)
*(S::RGBSpectrum, x::Number) = x*S
/(S::RGBSpectrum, x::Number) = RGBSpectrum(S.r/x, S.g/x, S.b/x)
*(a::RGBSpectrum, b::RGBSpectrum) = RGBSpectrum(a.r*b.r, a.g*b.g, a.b*b.b)
broadcast(f, S::RGBSpectrum) = RGBSpectrum(f(S.r), f(S.g), f(S.b))



# ------------------------------------------------
# Single line spectra

struct SingleLine <: Spectrum
    wavelength::Float64
    intensity::Float64
end

isblack(S::SingleLine) = S.intensity ≈ 0.0

# Addition of SingleLine objects
function +(S1::SingleLine, S2::SingleLine)
    if S1.wavelength != S2.wavelength
        error("Spectrum wavelengths don't match:\n$S1 \n $S2")
    end
    return SingleLine(S1.wavelength, S1.intensity + S2.intensity)
end
+(S::SingleLine, c::Real) = SingleLine(S.wavelength, S.intensity + c)
+(c::Real, S::SingleLine) = S+c
-(c::Real, S::SingleLine) = SingleLine(S.wavelength, S.intensity - c)


# Multiplication of SingleLine objects
function *(S1::SingleLine, S2::SingleLine)
    if S1.wavelength != S2.wavelength
        SingleLine(S1.wavelength, 0.0)
    end
    return SingleLine(S1.wavelength, S1.intensity * S2.intensity)
end

function *(S1::SingleLine, S2::SampledSpectrum)
    return SingleLine(S1.wavelength, S1.intensity * interpolate(S2, S1.wavelength))
end
*(S1::SampledSpectrum, S2::SingleLine) = S2*S1

*(c::Real, S::SingleLine) = SingleLine(S.wavelength, c*S.intensity)
*(S::SingleLine, c::Real) = c*S

/(S::SingleLine, c::Real) = SingleLine(S.wavelength, S.intensity/c)

interpolate(S::SingleLine, wavelength::Real) = S.wavelength ≈ wavelength ? S.intensity : 0.0


function to_XYZ(S::SingleLine)
    x,y,z = 0.0, 0.0, 0.0
    for i = 1:N_CIE-1
        if CIE_LAMBDA[i] < S.wavelength && CIE_LAMBDA[i+1] >= S.wavelength
            x += S.intensity * (CIE_X[i] + CIE_X[i+1]) / 2
            y += S.intensity * (CIE_Y[i] + CIE_Y[i+1]) / 2
            z += S.intensity * (CIE_Z[i] + CIE_Z[i+1]) / 2
            break
        end
    end
    return x, y, z
end

broadcast(f, S::SingleLine, y...) = SingleLine(S.wavelength, f(S.intensity, y...))



# ------------------------------------------------
# Conversion between XYZ and RGB

function XYZtoRGB(xyz::Array{Float64,1})
    return [3.240479*xyz[1] - 1.537150*xyz[2] - 0.498535*xyz[3],
           -0.969256*xyz[1] + 1.875991*xyz[2] + 0.041556*xyz[3],
            0.055648*xyz[1] - 0.204043*xyz[2] + 1.057311*xyz[3]]
end  

function RGBtoXYZ(rgb::Union{Array{Float64,1}, RGBSpectrum})
    return [0.412453*rgb[1] + 0.357580*rgb[2] + 0.180423*rgb[3],
            0.212671*rgb[1] + 0.715160*rgb[2] + 0.072169*rgb[3],
            0.019334*rgb[1] + 0.119193*rgb[2] + 0.950227*rgb[3]]
end
