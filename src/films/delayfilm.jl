

struct DelayFilm <: Film
    wavelength::Float64
    resX::Int64
    resY::Int64
    tres::Int64
    tmin::Float64
    tmax::Float64
    histogram::Array{Float64,1}
    counts::Array{Int64,1}
end

function DelayFilm(lmd::Real, resX::Integer, resY::Integer, tres::Integer, tmin::Real, tmax::Real)
    DelayFilm(lmd, resX, resY, tres, tmin, tmax, zeros(Float64, tres), zeros(Int64,tres))
end

uses_isect(F::DelayFilm) = true

function find_bin(t::Real, F::DelayFilm)
    nt = (t - F.tmin) / (F.tmax - F.tmin)
    return ceil(Int64, nt * F.tres)
end

function add_sample!(F::DelayFilm, sample::Sample, L::Spectrum, isect::Nullable{Intersection{T}}) where T<:Primitive
    if isblack(L) || isnull(isect)
        return nothing
    end
    isect = get(isect)
    if isect.tmin > F.tmax || isect.tmin < F.tmin
        return nothing
    end
    idx = find_bin(isect.tmin, F)
    I = interpolate(L, F.wavelength)
    F.histogram[idx] += I
    F.counts[idx] += 1
    return nothing
end

function write_txt(F::DelayFilm, fname::String)
    f = open(fname, "w")
    write(f, "# Distance and delay\n")
    write(f, "# lambda = $(F.wavelength)\n")
    write(f, "# tmin = $(F.tmin)\n")
    write(f, "# tmax = $(F.tmax)\n")
    write(f, "# tres = $(F.tres)\n")
    write(f, "# step = $((F.tmax - F.tmin) / F.tres)\n")
    for i = 1:F.tres
        value = F.counts[i] > 0 ? F.histogram[i] / F.counts[i] : 0.0
        write(f, string(value))
        write(f, "\n")
    end
    close(f)
end
