
"""
    DelayFilm

The delay film stores a 1D-histogram of observed intensity as a function of depth in the
scene. This can then be post-processed into e.g. reflectivity as a function of time delay.

"""
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


"""
    find_bin(F::DelayFilm, d::Real)

For a given distance `d`, find the corresponding bin in the `DelayFilm` histogram.

"""
function find_bin(F::DelayFilm, d::Real)
    nd = (d - F.tmin) / (F.tmax - F.tmin)
    return ceil(Int64, nd * F.nbins)
end



"""
    add_sample!(F::DelayFilm, sample::Sample, L::Spectrum, isect::Union{Intersection,Nothing})

Add a ray sample to the `DelayFilm`. The ray intersection `isect` is used to find the depth
in the scene where the reflection happened. The spectral intensity `L` will be added to the
bin corresponding to that distance. The camera sample `sample` is not used, but it's part
of the `Film` API.

"""
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


"""
    reset!(F::DelayFilm)

Reset all the pixel values of the `DelayFilm` to zero. This is mainly useful if you want to
use the same film in multiple computations, writing out the data and resetting in between.

"""
function reset!(F::DelayFilm)
    F.counts .= 0
    F.histogram .= 0.0
    nothing
end


"""
    write_txt(F::DelayFilm, fname::String)

Write the `DelayFilm` histogram into a text file. The file starts with seven lines of
metadata, followed by the values in the histogram vector.

# Metadata header example

```
# Distance and delay
# lambda = 532.0
# tmin = 0.0
# tmax = 4.0
# nbins = 1000
# step = 0.004
# counts = 2359296
```

"""
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
