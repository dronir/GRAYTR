
mutable struct PhotometricFilm{S<:Spectrum} <: Film
    integral::S
    N::Int64
    area::Float64
    resX::Int64
    resY::Int64
end

PhotometricFilm() = PhotometricFilm(nolight, 0)
PhotometricFilm(res::Integer) = PhotometricFilm(nolight, 0, res, res)

uses_isect(F::PhotometricFilm) = false


function add_sample!(F::PhotometricFilm, sample::Sample, L::Spectrum)
    F.N += 1
    isblack(L) && return nothing
    F.integral += L
    return nothing
end

function write_txt(F::PhotometricFilm, fname::String)
    S = F.area / F.N * F.integral
    
    return nothing
end
