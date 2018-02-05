
abstract type Spectrum end

struct SampledSpectrum
    low::Float64
    high::Float64
    values::Array{Float64,1}
end

import Base.+
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

