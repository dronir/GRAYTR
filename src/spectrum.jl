import Base.+, Base.-, Base.*, Base./, Base.convert, Base.zero, Base.one

export SampledSpectrum, SingleLine, NoLight, nolight

using Pkg
using DelimitedFiles


struct SpectrumStyle <: Broadcast.BroadcastStyle end
struct LineStyle <: Broadcast.BroadcastStyle end


# ------------------------------------------------
# Load CIE spectra


const CIE_PATH = joinpath(dirname(pathof(GRAYTR)), "..", "data", "lin2012xyz2e_1_7sf.csv")

const raw_CIE_data = readdlm(CIE_PATH, ',')
const N_CIE = size(raw_CIE_data)[1]

const CIE_DATA = SMatrix{4, N_CIE, Float64}(raw_CIE_data')
const CIE_YINT = sum(CIE_DATA[3,:])

#const CIE_LAMBDA = vec(raw_CIE_data[:,1])
#const CIE_X = vec(raw_CIE_data[:,2])
#const CIE_Y = vec(raw_CIE_data[:,3])
#const CIE_Z = vec(raw_CIE_data[:,4])

function make_CIE_functions(low, high, N)
end




# ------------------------------------------------
# The NoLight dummy spectrum

struct NoLight <: Spectrum end

const nolight = NoLight()

@inline isblack(N::NoLight) = true
@inline to_XYZ(S::NoLight) = (0.0, 0.0, 0.0)

Base.broadcastable(N::NoLight) = N

@inline interpolate(N::NoLight, x) = nolight
@inline integrate(N::NoLight) = 0.0

@inline Base.broadcasted(::typeof(+), s::NoLight, t::NoLight) = NoLight()
@inline Base.broadcasted(::typeof(*), n::NoLight, t::NoLight) = NoLight()

@inline Base.broadcasted(::typeof(*), n::Number, t::NoLight) = NoLight()
@inline Base.broadcasted(::typeof(*), t::NoLight, n::Number) = NoLight()
@inline Base.broadcasted(::typeof(+), n::Number, t::NoLight) = n
@inline Base.broadcasted(::typeof(+), t::NoLight, n::Number) = n

@inline Base.broadcasted(::typeof(*), t::NoLight, s::Spectrum) = NoLight()
@inline Base.broadcasted(::typeof(*), s::Spectrum, t::NoLight) = NoLight()
@inline Base.broadcasted(::typeof(+), t::NoLight, s::Spectrum) = s
@inline Base.broadcasted(::typeof(+), s::Spectrum, t::NoLight) = s

@inline Base.broadcasted(::typeof(*), t::NoLight, s::Broadcast.Broadcasted{SpectrumStyle,N}) where N = NoLight()
@inline Base.broadcasted(::typeof(+), t::NoLight, s::Broadcast.Broadcasted{SpectrumStyle,N}) where N = s




# ------------------------------------------------
# The SampledSpectrum type

struct SampledSpectrum <: Spectrum
    low::Int64
    high::Int64
    values::Array{Float64,1}
    delta::Float64
    invdelta::Float64
end


"""
    SampledSpectrum(low::Int64, high::Int64, values::Vector{Float64}) 

Create a new `SampledSpectrum` object with a wavelength range from `low` to `high`
nanometres and binned spectrum values given by `values`.

"""
function SampledSpectrum(low::Int64, high::Int64, values::Vector{Float64}) 
    delta = (high-low) / (length(values) - 1)
    SampledSpectrum(low, high, values, delta, 1/delta)
end


Base.maximum(S::SampledSpectrum) = Base.maximum(S.values)
Base.minimum(S::SampledSpectrum) = Base.minimum(S.values)

isblack(S::SampledSpectrum) = Base.maximum(S.values) ≈ 0.0



# Various functions
Base.length(S::SampledSpectrum) = Base.length(S.values)
Base.size(S::SampledSpectrum) = Base.size(S.values)
Base.getindex(S::SampledSpectrum, inds::Vararg{Int,N}) where N = S.values[inds...]
Base.getindex(S::SampledSpectrum, inds::CartesianIndex{1}) = S.values[inds]
Base.setindex!(S::SampledSpectrum, val, inds::Vararg{Int,N}) where N = S.values[inds...] = val
Base.setindex!(S::SampledSpectrum, val, inds::CartesianIndex{1}) = S.values[inds] = val

Base.collect(S::SampledSpectrum) = S.values

Base.broadcastable(S::SampledSpectrum) = S

# Broadcasting style for SampledSpectrum
Base.BroadcastStyle(::Type{<:SampledSpectrum}) = SpectrumStyle()
Base.BroadcastStyle(::SpectrumStyle, ::SpectrumStyle) = SpectrumStyle()
Base.BroadcastStyle(::SpectrumStyle, ::Broadcast.AbstractArrayStyle{0}) = SpectrumStyle()


function Base.similar(bc::Broadcast.Broadcasted{SpectrumStyle}, ::Type{ElType}) where ElType
#    flat = Broadcast.flatten(bc)
#    if !check_limits(flat.args)
#        error("Spectrum bounds mismatch")
#    end
    first = find_first(bc)
    return SampledSpectrum(first.low, first.high, similar(Array{ElType}, axes(bc)))
end

find_first(bc::Base.Broadcast.Broadcasted) = find_first(bc.args)
find_first(args::Tuple) = find_first(find_first(args[1]), Base.tail(args))
find_first(x) = x
find_first(a::SampledSpectrum, rest) = a
find_first(::Any, rest) = find_first(rest)

check_limits(x) = true
check_limits(args::Tuple{Any}) = true
check_limits(args::Tuple) = check_limits(args[1], check_limits(Base.tail(args)))
check_limits(S1::SampledSpectrum, S2::SampledSpectrum) = (S1.low == S2.low) && (S1.high == S2.high)
check_limits(S::SampledSpectrum, ::Any) = true
check_limits(::Any, ::Any) = true


function integrate(S::SampledSpectrum)
    dl = (S.high - S.low) / length(S.values)
    return dl * sum(S.values)
end

@inline function wavelength(S::SampledSpectrum, i::Integer)
    return S.low + i * S.delta
end

function interpolate(S::SampledSpectrum, a::Real)
    if a < S.low
        return S.values[1]
    elseif a >= S.high
        return S.values[end]
    end
    
    i = Int64(fld(a - S.low, S.delta))
    t = (a - i * S.delta) * S.invdelta
    return @inbounds lerp(t, S.values[i+1], S.values[i+2])
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
    x, y, z = 0.0, 0.0, 0.0
    for i = 1:N_CIE
        lam = CIE_DATA[1,i]
        s = interpolate(S, lam)
        x += s * CIE_DATA[2,i]
        y += s * CIE_DATA[3,i]
        z += s * CIE_DATA[4,i]
    end
    return x/CIE_YINT, y/CIE_YINT, z/CIE_YINT
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

"""
    SingleLine(wavelength, value)

Create a `SingleLine` spectrum representing a delta-distribution of a value at a certain
wavelength.

"""
struct SingleLine <: Spectrum
    wavelength::Float64
    value::Float64
end

isblack(S::SingleLine) = S.value ≈ 0.0

# Addition of SingleLine objects
function +(S1::SingleLine, S2::SingleLine)
    if S1.wavelength != S2.wavelength
        error("Spectrum wavelengths don't match:\n$S1 \n $S2")
    end
    return SingleLine(S1.wavelength, S1.value + S2.value)
end
+(S::SingleLine, c::Real) = SingleLine(S.wavelength, S.value + c)
+(c::Real, S::SingleLine) = S+c
-(c::Real, S::SingleLine) = SingleLine(S.wavelength, S.value - c)


# Multiplication of SingleLine objects
function *(S1::SingleLine, S2::SingleLine)
    if S1.wavelength != S2.wavelength
        SingleLine(S1.wavelength, 0.0)
    end
    return SingleLine(S1.wavelength, S1.value * S2.value)
end


Base.broadcastable(S::SingleLine) = S

# Multiplication broadcasting

@inline function Base.broadcasted(::typeof(*), S1::SingleLine, S2::SingleLine)
    if S1.wavelength != S2.wavelength
        return SingleLine(S1.wavelength, 0.0)
    end
    return SingleLine(S1.wavelength, S1.value * S2.value)
end

@inline function Base.broadcasted(::typeof(*), S1::SingleLine, S2::SampledSpectrum)
    return SingleLine(S1.wavelength, S1.value * interpolate(S2, S1.wavelength))
end
@inline function Base.broadcasted(::typeof(*), S2::SampledSpectrum, S1::SingleLine)
    return SingleLine(S1.wavelength, S1.value * interpolate(S2, S1.wavelength))
end
@inline function Base.broadcasted(::typeof(*), S1::SingleLine, n::Number)
    return SingleLine(S1.wavelength, S1.value * n)
end
@inline function Base.broadcasted(::typeof(*), n::Number, S1::SingleLine)
    return SingleLine(S1.wavelength, S1.value * n)
end

@inline function Base.broadcasted(::typeof(/), S1::SingleLine, n::Number)
    return SingleLine(S1.wavelength, S1.value / n)
end

@inline function Base.broadcasted(::typeof(+), S1::SingleLine, S2::SingleLine)
    if S1.wavelength != S2.wavelength
        return SingleLine(S1.wavelength, 0.0)
    end
    return SingleLine(S1.wavelength, S1.value + S2.value)
end

@inline function Base.broadcasted(::typeof(+), S1::SingleLine, S2::SampledSpectrum)
    return SingleLine(S1.wavelength, S1.value + interpolate(S2, S1.wavelength))
end
@inline function Base.broadcasted(::typeof(+), S2::SampledSpectrum, S1::SingleLine)
    return SingleLine(S1.wavelength, S1.value + interpolate(S2, S1.wavelength))
end
@inline function Base.broadcasted(::typeof(+), S1::SingleLine, n::Number)
    return SingleLine(S1.wavelength, S1.value + n)
end
@inline function Base.broadcasted(::typeof(+), n::Number, S1::SingleLine)
    return SingleLine(S1.wavelength, S1.value + n)
end


@inline Base.broadcasted(::typeof(*), S::SingleLine, B::Broadcast.Broadcasted{SpectrumStyle,N}) where N = S .* materialize(B)
@inline Base.broadcasted(::typeof(+), S::SingleLine, B::Broadcast.Broadcasted{SpectrumStyle,N}) where N = S .+ materialize(B)




function *(S1::SingleLine, S2::SampledSpectrum)
    return SingleLine(S1.wavelength, S1.value * interpolate(S2, S1.wavelength))
end
*(S1::SampledSpectrum, S2::SingleLine) = S2*S1

*(c::Real, S::SingleLine) = SingleLine(S.wavelength, c*S.value)
*(S::SingleLine, c::Real) = c*S

/(S::SingleLine, c::Real) = SingleLine(S.wavelength, S.value/c)

interpolate(S::SingleLine, wavelength::Real) = S.wavelength ≈ wavelength ? S.value : 0.0

integrate(S::SingleLine) = S.value

function to_XYZ(S::SingleLine)
    x,y,z = 0.0, 0.0, 0.0
    for i = 1:N_CIE-1
        if CIE_DATA[1,i] < S.wavelength && CIE_DATA[1,i+1] >= S.wavelength
            t = (CIE_DATA[1,i+1] - S.wavelength) / (CIE_DATA[1,i+1] - CIE_DATA[1,i])
            x += S.value * lerp(t, CIE_DATA[2,i], CIE_DATA[2,i+1])
            y += S.value * lerp(t, CIE_DATA[3,i], CIE_DATA[3,i+1])
            z += S.value * lerp(t, CIE_DATA[4,i], CIE_DATA[4,i+1])
            break
        end
    end
    return x/CIE_YINT, y/CIE_YINT, z/CIE_YINT
end

broadcast(f, S::SingleLine, y...) = SingleLine(S.wavelength, f(S.value, y...))



# ------------------------------------------------
# Conversion between XYZ and RGB

function XYZtoRGB(x, y, z)
    return [3.240479*x - 1.537150*y - 0.498535*z,
           -0.969256*x + 1.875991*y + 0.041556*z,
            0.055648*x - 0.204043*y + 1.057311*z]
end  

function RGBtoXYZ(rgb::Union{Array{Float64,1}, RGBSpectrum})
    return [0.412453*rgb[1] + 0.357580*rgb[2] + 0.180423*rgb[3],
            0.212671*rgb[1] + 0.715160*rgb[2] + 0.072169*rgb[3],
            0.019334*rgb[1] + 0.119193*rgb[2] + 0.950227*rgb[3]]
end
