using FileIO
using Images




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



"""
    ImageFilm

A film that produces a resolved image. The resolution of the image is `(resX, resY)`,
the `gain` value will be used 
"""
struct ImageFilm{T<:Filter} <: Film
    resX::Int64
    resY::Int64
    gain::Float64
    filter::T
    pixels::Array{Pixel,2}
    filtertable::Array{Float64,2}
    CIE_table::CIE_Table
end




const FILTERTABLE_SIZE = 16

function make_filtertable(f::Filter)
    ftbl = zeros(Float64, (FILTERTABLE_SIZE, FILTERTABLE_SIZE))
    for j = 1:FILTERTABLE_SIZE
        fy = (j + 0.5) * f.ywidth / FILTERTABLE_SIZE
        for i = 1:FILTERTABLE_SIZE
            fx = (i + 0.5) * f.xwidth / FILTERTABLE_SIZE
            ftbl[i,j] = evaluate(f, fx, fy)
        end
    end
    ftbl
end

ImageFilm(x::Integer, y::Integer, gain::Real, f::Filter, CIE_table::CIE_Table) = ImageFilm(x, y, gain, f, 
                                                            zeros(Pixel,(x,y)), 
                                                            make_filtertable(f), CIE_table::CIE_Table)


"""
    add_sample!(F::ImageFilm, sample::Sample, L::Spectrum, isect::Union{Intersection,Nothing})

Add a ray sample to the `ImageFilm`. The image coordinates in `sample` tell which pixel we
are adding to. The spectral intensity `L` will be added to the pixel, and to nearby pixels
depending on the filter included in the film. `isect` is not used, but it's part of the
`Film` API.

"""
function add_sample!(F::ImageFilm, sample::Sample, L::Spectrum, isect::Union{Intersection,Nothing})
    dimgX = sample.imgX - 0.5
    dimgY = sample.imgY - 0.5
    x0 = ceil(Int64, dimgX - F.filter.xwidth)
    x1 = floor(Int64, dimgX + F.filter.xwidth)
    y0 = ceil(Int64, dimgY - F.filter.ywidth)
    y1 = floor(Int64, dimgY + F.filter.ywidth)
    x0 = max(x0, 1)
    x1 = min(x1, F.resX)
    y0 = max(y0, 1)
    y1 = min(y1, F.resY)
    if (x1-x0) < 0 || (y1-y0) < 0
        return
    end
    x, y, z = to_XYZ(L, F.CIE_table)
    for i = x0:x1
        for j = y0:y1
            # TODO: unoptimized way, without using filter table
            w = evaluate(F.filter, (i-dimgX), (j-dimgY))
            @inbounds F.pixels[i,j] = add(F.pixels[i,j], w*x, w*y, w*z, w)
        end
    end
end


"""
    reset!(F::ImageFilm)

Reset all the pixel values of the `ImageFilm` to zero. This is mainly useful if you want to
use the same `ImageFilm` in multiple computations, writing out the image and resetting in
between.

"""
function reset!(F::ImageFilm)
    for i = 1:F.resX
        for j = 1:F.resY
            F.pixels[i,j] = zero(Pixel)
        end
    end
    nothing
end


"""
    write_image(F::ImageFilm, fname::String="test.png")

Write the given image into a file determined by the filename.
"""
function write_image(F::ImageFilm, fname::String="test.png")
    data = zeros(Float64, (3, size(F.pixels)...))
    for i = 1:size(F.pixels,1)
        for j = 1:size(F.pixels,2)
            P = F.pixels[i,j]
            r, g, b = XYZtoRGB(F.gain*P.x/P.w, F.gain*P.y/P.w, F.gain*P.z/P.w) 
            data[1,i,j] = r
            data[2,i,j] = g
            data[3,i,j] = b
            
        end
    end
    img = colorview(RGB, map(clamp01nan, data))
    save(fname, img)
end

