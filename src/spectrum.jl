import Base.+, Base.-, Base.*, Base./, Base.convert, Base.zero, Base.one

export SampledSpectrum, SingleLine, NoLight, nolight, make_CIE_table

using Pkg
using DelimitedFiles


struct SpectrumStyle <: Broadcast.AbstractArrayStyle{1} end
struct LineStyle <: Broadcast.BroadcastStyle end


# ------------------------------------------------
# Load CIE spectra


const CIE_PATH = joinpath(dirname(pathof(GRAYTR)), "..", "data", "lin2012xyz2e_1_7sf.csv")

const raw_CIE_data = readdlm(CIE_PATH, ',')
const N_CIE = size(raw_CIE_data)[1]

const CIE_DATA = SMatrix{4, N_CIE, Float64}(raw_CIE_data')
const CIE_YINT = sum(CIE_DATA[3,:])


struct CIE_Table
    data::Array{Float64,2}
    yint::Float64
end

CIE_Table(data) = CIE_Table(data, sum(data[:,2]))




# ------------------------------------------------
# Utility functions used when integrating spectral data

function find_bin(data, x)
    N_DATA = size(data)[1]
    for i = 1:N_DATA-1
        if data[i,1] <= x && data[i+1,1] > x
            return i
        end
    end
    nothing
end


function get_midpoint(DATA, x, n)
    a = DATA[n,1]
    b = DATA[n+1,1]
    t = (x - a) / (b - a)
    y = lerp(t, DATA[n,2], DATA[n+1,2])
    return t,y
end


function integrate(X, Y)
    deltas = X[2:end] - X[1:end-1]
    means = (Y[1:end-1] + Y[2:end]) / 2
    return sum(deltas .* means)
end


"""
    compute_bin(DATA; low, high ; averaged = false)

Compute the bin value by integrating over a (hopefully) higher-resolution spectrum.

- `DATA` is a `(N, 2)` array with wavelengths (in nm) in the first column and values in the second.
- `low` and `high` are the wavelength limits of the bin.
- `averaged` optionally returns an average over the spectrum instead of an integral 
  (i.e. it divides the result by the bin length).

"""
function compute_bin(DATA, low, high ; averaged = false)
    nlow = find_bin(DATA, low)
    nhigh = find_bin(DATA, high)
    
    if nlow == nothing || nhigh == nothing
        return 0.0
    end
    
    if nlow == nhigh
        # There is only one data value in the range.
        # Return it if average, or return it times the interval if non-averaged
        low_t, low_y = get_midpoint(DATA, low, nlow)
        high_t, high_y = get_midpoint(DATA, high, nhigh)
        return 0.5 * (low_y + high_y) * (high - low)
    end
    
    low_t, low_y = get_midpoint(DATA, low, nlow)
    high_t, high_y = get_midpoint(DATA, high, nhigh)
    
    
    low_part =  0.5 * (low_y + DATA[nlow+1,2]) * (DATA[nlow+1,1] - low)
    high_part = 0.5 * (high_y + DATA[nlow,2]) * (high - DATA[nhigh,1])
    
    mid_part = integrate(DATA[nlow+1 : nhigh, 1], DATA[nlow+1 : nhigh, 2])

    S = low_part + mid_part + high_part

    return averaged ? S / (high - low) : S
end





# ------------------------------------------------
# The NoLight dummy spectrum

struct NoLight <: Spectrum end

const nolight = NoLight()

@inline isblack(N::NoLight) = true
@inline to_XYZ(S::NoLight) = (0.0, 0.0, 0.0)
to_XYZ(S::NoLight, C::CIE_Table) = (0.0, 0.0, 0.0)

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
end


"""
    SampledSpectrum(low::Int64, high::Int64, values::Vector{Float64}) 

Create a new `SampledSpectrum` object with a wavelength range from `low` to `high`
nanometres and binned spectrum values given by `values`.

"""
function SampledSpectrum(low::Int64, high::Int64, values::Vector{Float64}) 
    delta = (high-low) / length(values)
    SampledSpectrum(low, high, values, delta)
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
Base.BroadcastStyle(::SpectrumStyle, ::Broadcast.DefaultArrayStyle{0}) = SpectrumStyle()


function Base.similar(bc::Broadcast.Broadcasted{SpectrumStyle}, ::Type{ElType}) where ElType
    flat = Broadcast.flatten(bc)
    first = find_first(flat)
    l = first.low
    h = first.high
    N = length(first.values)
    
    for x in flat.args
        if isa(x, SampledSpectrum)
            if (x.low != l) || (x.high != h) || (length(x.values) != N)
                error("Spectral shape mismatch!")
            end
        end
    end
    
    return SampledSpectrum(first.low, first.high, similar(Array{ElType}, axes(bc)))
end

find_first(bc::Base.Broadcast.Broadcasted) = find_first(bc.args)
find_first(args::Tuple) = find_first(find_first(args[1]), Base.tail(args))
find_first(args::Tuple{}) = nothing
find_first(x) = x
find_first(a::SampledSpectrum, rest) = a
find_first(::Any, rest) = find_first(rest)


@inline isSS(S::Any) = false
@inline isSS(S::SampledSpectrum) = true



integrate(S::SampledSpectrum) = sum(S.values)


function interpolate(S::SampledSpectrum, a::Real)
    if a <= S.low || a >= S.high
        return 0.0
    end
    
    i = Int64(fld(a - S.low, S.delta))
    return S.values[i+1]
#    t = (a - i * S.delta) / S.delta
#    return @inbounds lerp(t, S.values[i+1], S.values[i+2])
end



function make_CIE_table(S::SampledSpectrum)
    N = length(S)
    data = zeros(N,3)
    delta = (S.high - S.low) / N
    CIE_lambda = raw_CIE_data[:,1]
    for n = 1:N
        # (low_x, high_x) are the wavelength limits of the bin; compute average over that
        low_x = S.low + (n-1) * delta
        high_x = S.low + n * delta
        data[n,1] = compute_bin(hcat(CIE_lambda, raw_CIE_data[:,2]), low_x, high_x ; averaged = true)
        data[n,2] = compute_bin(hcat(CIE_lambda, raw_CIE_data[:,3]), low_x, high_x ; averaged = true)
        data[n,3] = compute_bin(hcat(CIE_lambda, raw_CIE_data[:,4]), low_x, high_x ; averaged = true)
    end
    return CIE_Table(data)
end



function to_XYZ(S::SampledSpectrum, C::CIE_Table)
    x =y = z = 0.0
    for i = 1:length(S)
        x += S.values[i] * C.data[i,1]
        y += S.values[i] * C.data[i,2]
        z += S.values[i] * C.data[i,3]
    end
    return x/C.yint, y/C.yint, z/C.yint
end





function to_XYZ(S::SampledSpectrum)
    x, y, z = 0.0, 0.0, 0.0
    @inbounds for i = 1:N_CIE
        lam = CIE_DATA[1,i]
        s = interpolate(S, lam)
        x += s * CIE_DATA[2,i]
        y += s * CIE_DATA[3,i]
        z += s * CIE_DATA[4,i]
    end
    return x/CIE_YINT, y/CIE_YINT, z/CIE_YINT
end






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


const Adobe_RGB = reshape([2.0413690, -0.5649464, -0.3446944,-0.9692660,  1.8760108,  0.0415560, 0.0134474, -0.1183897,  1.0154096], (3,3))'
const CIE_RGB = reshape([2.3706743, -0.9000405, -0.4706338,-0.5138850,  1.4253036,  0.0885814, 0.0052982, -0.0146949,  1.0093968], (3,3))'


# ------------------------------------------------
# Conversion between XYZ and RGB

XYZtoRGB(xyz::Vector) = CIE_RGB * xyz


function XYZtoRGB(x, y, z)
    return 2.3706743*x - 0.9000405*y - 0.4706338*z,
           -0.5138850*x + 1.4253036*y + 0.0885814*z,
            0.0052982*x - 0.0146949*y + 1.0093968*z
end  
