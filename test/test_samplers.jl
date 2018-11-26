using Base.Test
include("abstracts.jl")
include("samplers.jl")

struct DummySampler <: Sampler
    xstart::Int64
    xend::Int64
    ystart::Int64
    yend::Int64
end

@testset "Subwindow computation" begin
    D = DummySampler(1, 32, 1, 32)
    @test compute_subwindow(D, 1, 4) == (1, 16, 1, 16)
    @test compute_subwindow(D, 2, 4) == (17, 32, 1, 16)
    @test compute_subwindow(D, 3, 4) == (1, 16, 17, 32)
    @test compute_subwindow(D, 4, 4) == (17, 32, 17, 32)
end

nothing
