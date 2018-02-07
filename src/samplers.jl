

# Sampler interface:
#  get_samples()
#  maximumsamplecount()
#  getsubsampler()

# TODO MOVE ELSEWHERE
lerp(t, a, b) = (1-t)*a + t*b


# Stratified1D and Stratified2D generate a bunch of stratified samples in the [0,1] interval.

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
    dx = sampler.xend - sampler.xstart + 1
    dy = sampler.yend - sampler.ystart + 1
    n -= 1
    nx = count
    ny = 1
    while nx%2==0 && 2*dx*ny < dy*nx
        nx = div(nx, 2)
        ny *= 2
    end

    cx = div(dx,nx)
    cy = div(dy,ny)

    x0 = n % nx 
    y0 = div(n, ny) 
    return x0*cx+1, (x0+1)*cx, y0*cy+1, (y0+1)*cy
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

StratifiedSampler(x0::Int64, x1::Int64, y0::Int64, y1::Int64, xs::Int64, ys::Int64, j::Bool) = StratifiedSampler(x0,x1,y0,y1,xs,ys,xs*ys,j)

roundsize(S::StratifiedSampler, size::Integer) = size


function get_subsampler(sampler::StratifiedSampler, n::Integer, count::Integer)
    x0, x1, y0, y1 = compute_subwindow(sampler, n, count)
    if x0 == x1 || y0 == y1
        return Nullable{StratifiedSampler}()
    else
        return Nullable(StratifiedSampler(x0, x1, y0, y1, sampler.xs, sampler.ys, 
                                          sampler.jitter))
    end
end

# Returns an empty array when all the sampling has been done
function get_samples(sampler::StratifiedSampler, state::Integer)
    dx = sampler.xend - sampler.xstart + 1
    dy = sampler.yend - sampler.ystart + 1
    if state >= dx*dy
        return CameraSample[], state
    end

    xpos = sampler.xstart + state % dx
    ypos = sampler.ystart + div(state, dy)
    N = sampler.xs * sampler.ys
    out = Array{CameraSample}(N)
    img_samples = stratified2D(N, sampler.jitter)
    lens_samples = stratified2D(N, sampler.jitter)
    shuffle!(lens_samples[:,1])
    shuffle!(lens_samples[:,2])
    for i = 1:N
        s = CameraSample(
            img_samples[i,1] + xpos,
            img_samples[i,2] + ypos,
            lens_samples[i,1],
            lens_samples[i,2]
        )
        out[i] = s
    end
    return out, state+1
end
