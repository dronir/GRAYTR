
"""
    PhotometricFilm

The photometric film is a one-pixel detector, simply adding together all the intensity
coming from the scene to the camera.

TODO: This is really not in working condition, do not use it.

"""
mutable struct PhotometricFilm{S<:Spectrum} <: Film
    integral::S
    N::Int64
    area::Float64
    resX::Int64
    resY::Int64
end


"""
    add_sample!(F::PhotometricFilm, sample::Sample, L::Spectrum, isect::Union{Intersection,Nothing})

Add a ray sample to the `ImageFilm`. The spectral intensity `L` will be added to the
detector's single pixel. The ray intersection `isect` nor the camera sample `sample` are
not used, but they are part of the `Film` API.

TODO: This relies heavily on the addition of Spectrum types, which is a bit poorly defined
currently.

"""
function add_sample!(F::PhotometricFilm, sample::Sample, L::Spectrum, isect::Union{Intersection,Nothing})
    F.N += 1
    isblack(L) && return nothing
    F.integral += L
    return nothing
end


"""
    write_txt(F::PhotometricFilm, fname::String)

Write the value of the single pixel of the `PhotometricFilm` into a file named `fname`.

"""
function write_txt(F::PhotometricFilm, fname::String)
    S = F.area / F.N * F.integral    
    f = open(fname, "w")
    write(f, string(S))
    return nothing
end
