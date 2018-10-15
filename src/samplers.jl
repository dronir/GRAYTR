

"""
    stratified1D(N::Integer, jitter::Bool)

Divides the `[0,1]` interval into `N` bins and either returns an array with the bin centers 
(if `jitter == false`), or a uniformly random point inside each bin (if `jitter == true`).

# Examples
```julia-repl
julia> stratified1D(2, false)
2-element Array{Float64,1}:
 0.25
 0.75
```
"""
stratified1D(N::Integer, jitter::Bool) = [(i-1 + (jitter ? rand() : 0.5))/N for i = 1:N]


"""
    stratified2D(N::Integer, jitter::Bool)

Returns a two-column array with each column containing values of the [0,1] interval divided
into `N` bins, with the values at the bin centres if `jitter == false` and randomly
jittered inside the bin if `jitter == true`.

"""
function stratified2D(N::Integer, jitter::Bool)
    out = zeros(Float64, (N,2))
    for i = 1:N
        out[i,1] = (i-1 + (jitter ? rand() : 0.5))/N
        out[i,2] = (i-1 + (jitter ? rand() : 0.5))/N
    end
    out
end


"""
    compute_subwindow(sampler::Sampler, n::Int64, count::Int64)

This computes the pixel range which a given subsampler samples, given the index `n` of the
subsampler and the total number of areas, `count`.

Returns four integers, given the pixel coordinates (x0, x1, y0, y1) of the subwindow
(inclusive).

Note that indexing starts from 1 unlike in book's C++ code.

"""
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


"""
    CameraSample

A wrapper for four floating point values, generally to represent two random points needed
for ray creation from a camera: one point on the film and the other on the lens.

"""
struct CameraSample <: Sample
    imgX::Float64
    imgY::Float64
    lensU::Float64
    lensV::Float64
end




"""
    StratifiedSampler <: Sampler

A StratifiedSampler will generate CameraSamples. It handles a certain window of pixels and
can be divided into subsamplers which are new copies of StratifiedSampler with the desired
subwindow of pixels.

It can generate multiple samples per pixel, either divided into clean bins, or jittered.

# Fields

- `xstart::Int64`: The start of the pixel range in the x direction.
- `xend::Int64`: The end of the pixel range in the x direction.
- `ystart::Int64`: The start of the pixel range in the y direction.
- `yend::Int64`: The end of the pixel range in the y direction.
- `xs::Int64`: The number of sample bins per pixel in the x direction.
- `ys::Int64`: The number of sample bins per pixel in the y direction.
- `jitter::Bool`: Whether to jitter samples or not.

"""
struct StratifiedSampler <: Sampler
    xstart::Int64
    xend::Int64
    ystart::Int64
    yend::Int64
    xs::Int64
    ys::Int64
    jitter::Bool
end

StratifiedSampler(x0::Int64, x1::Int64, y0::Int64, y1::Int64, xs::Int64, ys::Int64, 
                  j::Bool) = StratifiedSampler(x0,x1,y0,y1,xs,ys,j)

StratifiedSampler(rangeX::Int64, rangeY::Int64, res::Int64) = StratifiedSampler(1,rangeX,1,
                                                               rangeY,res,res,true)


"""
    get_subsampler(sampler::StratifiedSampler, n::Integer, count::Integer)

Given a `StratifiedSampler`, a total number of subwindows `count` and a subwindow index
`n`, returns a new `StratifiedSampler` which samples that subwindow.

"""
function get_subsampler(sampler::StratifiedSampler, n::Integer, count::Integer)
    x0, x1, y0, y1 = compute_subwindow(sampler, n, count)
    if x0 == x1 || y0 == y1
        return nothing
    else
        return StratifiedSampler(x0, x1, y0, y1, sampler.xs, sampler.ys, sampler.jitter)
    end
end


"""
    finished(sampler::StratifiedSampler, state::Integer)

Returns `true` when `state` indicates that `sampler` is "finished", i.e. cannot return
any more samples.

"""
function finished(sampler::StratifiedSampler, state::Integer)
    dx = sampler.xend - sampler.xstart + 1
    dy = sampler.yend - sampler.ystart + 1
    return state >= dx*dy
end


"""
    get_samples(sampler::StratifiedSampler, state::Integer)

Generate samples from the `StratifiedSampler`, given its `state`, Returns the array of
samples and the new state of the sampler.

"""
function get_samples(sampler::StratifiedSampler, state::Integer)
    out = Array{CameraSample}(N)
    state = get_samples!(sampler, state, out)
    return out, state
end


"""
    get_samples!(sampler::StratifiedSampler, state::Integer, out::Array{CameraSample,1})

Generate samples for a single pixel from the `StratifiedSampler`, given its `state`, and
fill the array `out` with them. Returns the updated state of the sampler.

"""
function get_samples!(sampler::StratifiedSampler, state::Integer, out::Array{CameraSample,1})
    dx = sampler.xend - sampler.xstart + 1
    dy = sampler.yend - sampler.ystart + 1
    if state >= dx*dy
        return state
    end

    xpos = sampler.xstart + state % dx
    ypos = sampler.ystart + div(state, dy)
    N = sampler.xs * sampler.ys
    img_samples = stratified2D(N, sampler.jitter)
    lens_samples = stratified2D(N, sampler.jitter)
    shuffle!(lens_samples[:,1])
    shuffle!(lens_samples[:,2])
    for i = 1:N
        out[i] = CameraSample(
            img_samples[i,1] + xpos,
            img_samples[i,2] + ypos,
            lens_samples[i,1],
            lens_samples[i,2]
        )
    end
    return state+1
end
