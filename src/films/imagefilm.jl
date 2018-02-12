struct ImageFilm{T<:Filter} <: Film
    resX::Int64
    resY::Int64
    filter::T
    pixels::Array{Pixel,2}
    filtertable::Array{Float64,2}
end

uses_isect(F::ImageFilm) = false

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

ImageFilm(x::Integer, y::Integer, f::Filter) = ImageFilm(x, y, f, zeros(Pixel,(x,y)), 
                                                         make_filtertable(f))

function add_sample!(F::ImageFilm, sample::Sample, L::Spectrum)
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
    xyz = to_XYZ(L)
    for i = x0:x1
        for j = y0:y1
            # TODO: unoptimized way, without using filter table
            w = evaluate(F.filter, (i-dimgX), (j-dimgY))
            F.pixels[i,j] = add(F.pixels[i,j], w*xyz, w)
        end
    end
end

function write_image(F::ImageFilm) 
    data = zeros(Float64, (3, size(F.pixels)...))
    for i = 1:size(F.pixels,1)
        for j = 1:size(F.pixels,2)
            P = F.pixels[i,j]
            data[:,i,j] = XYZtoRGB([P.x/P.w, P.y/P.w, P.z/P.w])
        end
    end
    img = colorview(RGB, map(clamp01nan, data))
    save("test.png", img)
end

function write_bwtxt(F::ImageFilm) 
    f = open("test.txt", "w")
    for i = 1:minimum(size(F.pixels))
        write(f, repr(F.pixels[i,i]))
        write(f, "\n")
    end
end
