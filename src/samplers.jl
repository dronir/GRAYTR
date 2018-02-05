
abstract type Sampler end
abstract type Sample end

# Sampler interface:
#  get_samples()
#  maximumsamplecount()
#  getsubsampler()

# TODO CHECK CORRECTNESS
# TODO MOVE ELSEWHERE
lerp(t, a, b) = a + t * (b-a)

stratified1D(N::Integer, jitter::Bool) = [(i-1 + (jitter ? rand() : 0.5))/N for i = 0:N:1]
function stratified2D(N::Integer, jitter::Bool)
    out = zeros(Float64, (N,2))
    for i = 1:N
        out[i,1] = (i-1 + (jitter ? rand() : 0.5))/N
        out[i,2] = (i-1 + (jitter ? rand() : 0.5))/N
    end
    out
end

# This computes the pixel range which a given sampler samples, given
# the number `n` of the area and the total number `count` of areas.
# Note that indexing starts from 1 unlike in book's C++ code.
function compute_subwindow(sampler::Sampler, n::Int64, count::Int64)
    xres = sampler.xend - sampler.xstart
    yres = sampler.yend - sampler.ystart
    n -= 1
    nx = count
    ny = 1
    while nx%2==0 && 2*xres*ny < yres*nx
        nx = div(nx, 2)
        ny *= 2
    end
    x0 = n % nx + 1
    y0 = div(n, ny) + 1
    tx0 = x0 / nx
    ty0 = t0 / ny
    tx1 = (x0+1) / nx
    ty1 = (y0+1) / ny
    xstart1 = floor(Int64, lerp(tx0, sampler.xstart, sampler.xend))
    xendt1 = floor(Int64, lerp(tx0, sampler.xstart, sampler.xend))
    ystart1 = floor(Int64, lerp(tx0, sampler.ystart, sampler.yend))
    yendt1 = floor(Int64, lerp(tx0, sampler.ystart, sampler.yend))
end


struct CameraSample <: Sample
    imgX::Float64
    imgY::Float64
    lensU::Float64
    lensV::Float64
end

# (xstart, xend), (ystart, yend) is the pixel range handled by the sampler
# (xs, ys) is the number of strata in each direction
struct StratifiedSampler <: Sampler
    xstart::Int64
    xend::Int64
    ystart::Int64
    yend::Int64
    xs::Int64
    ys::Int64
    samplesperpixel::Int64
    jitter::Bool
end

roundsize(S::StratifiedSampler, size::Integer) = size


function get_subsampler(sampler::StratifiedSampler, n::Integer, count::Integer)
    x0, x1, y0, y1 = compute_subwindow(sampler, n, count)
    if x0 == x1 || y0 == y1
        return Nullable{StratifiedSampler}()
    else
        return Nullable(StratifiedSampler(x0, x1, t0, t1, sampler.xs, sampler.ys, 
                                          sampler.samplesperpixel, sampler.jitter))
    end
end

# Returns an empty array when all the sampling has been done
function get_samples(sampler::StratifiedSampler, state::Integer)
    dx = sampler.xend - sampler.xstart
    dy = sampler.yend - sampler.ystart
    if state >= dx*dy
        return CameraSample[], state
    end

    xpos = sampler.xstart + state % dx
    ypos = sampler.ystart + div(state, dy)
    N = sampler.xs * sampler.xy
    out = Array{CameraSample}(N)
    img_samples = stratified2D(N, sampler.jitter)
    lens_samples = stratified2D(N, sampler.jitter)
    shuffle!(lens_samples[:,1])
    shuffle!(lens_samples[:,2])
    for i = 1:N
        out[i].imgX = img_samples(i,1)
        out[i].imgY = img_samples(i,2)
        out[i].lensU = lens_samples(i,1)
        out[i].lensV = lens_samples(i,2)
    end
    return out, state+1
end
