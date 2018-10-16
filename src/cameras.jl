

include("cameras/projective.jl")

"""
    Pixel

Representation of a pixel. A wrapper for the XYZ colors and a weight.
"""
struct Pixel
    x::Float64
    y::Float64
    z::Float64
    w::Float64
end

import Base.zero
zero(::Type{Pixel}) = Pixel(0.0,0.0,0.0,0.0)
add(P::Pixel, S::Array{Float64,1}) = Pixel(P.x+S[1], P.y+S[2], P.z+S[3], P.w+S[4])
add(P::Pixel, S::Array{Float64,1}, w::Float64) = Pixel(P.x+S[1], P.y+S[2], P.z+S[3], P.w+w)
add(P::Pixel, x, y, z, w) = Pixel(P.x+x, P.y+y, P.z+z, P.w+w)


