

include("cameras/projective.jl")

struct Pixel
    x::Float64
    y::Float64
    z::Float64
    w::Float64
end

import Base.zero
zero(::Type{Pixel}) = Pixel(0.0,0.0,0.0,0.0)
add(P::Pixel, S::Array{Float64,1}, w::Float64) = Pixel(P.x+S[1], P.y+S[2], P.z+S[3], P.w+w)

include("films/imagefilm.jl")
