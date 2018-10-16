

struct DelayFilm <: Film
    wavelength::Float64
    area::Float64
    resX::Int64
    resY::Int64
    nbins::Int64
    tmin::Float64
    tmax::Float64
    histogram::Array{Float64,1}
    counts::Array{Int64,1}
end

function DelayFilm(lmd::Real, window::Array{Float64,1}, resX::Integer, resY::Integer, nbins::Integer, tmin::Real, tmax::Real)
    area = (window[2] - window[1]) * (window[4] - window[3])
    DelayFilm(lmd, area, resX, resY, nbins, tmin, tmax, zeros(Float64, nbins), zeros(Int64,nbins+1))
end

function find_bin(F::DelayFilm, t::Real)
    nt = (t - F.tmin) / (F.tmax - F.tmin)
    return ceil(Int64, nt * F.nbins)
end

function add_sample!(F::DelayFilm, sample::Sample, L::Spectrum, isect::Union{Intersection,Nothing})
    F.counts[end] += 1
    if isect == nothing
        return nothing
    end   
    if isect.tmin > F.tmax || isect.tmin < F.tmin
        return nothing
    end
    idx = find_bin(F, isect.tmin)
    F.counts[idx] += 1
    
    if isblack(L)
        return nothing
    end
    I = interpolate(L, F.wavelength)
    F.histogram[idx] += I
    return nothing
end

function write_txt(F::DelayFilm, fname::String)
    f = open(fname, "w")
    write(f, "# Distance and delay\n")
    write(f, "# lambda = $(F.wavelength)\n")
    write(f, "# tmin = $(F.tmin)\n")
    write(f, "# tmax = $(F.tmax)\n")
    write(f, "# nbins = $(F.nbins)\n")
    write(f, "# step = $((F.tmax - F.tmin) / F.nbins)\n")
    write(f, "# counts = $(F.counts[end])\n")
    sum_counts = F.counts[end]
    for i = 1:F.nbins
        area_in_bin = F.area * F.counts[i] / sum_counts
        value = F.counts[i] > 0 ? area_in_bin * F.histogram[i] / F.counts[i] : 0.0
        write(f, string(value))
        write(f, "\n")
    end
    close(f)
end
