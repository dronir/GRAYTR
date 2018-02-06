using Base.Test
include("abstracts.jl")
include("geometry.jl")
using Geometry
include("samplers.jl")
include("rays.jl")
include("cameras.jl")
include("filters.jl")

struct DummyFilm <: Film 
    xres::Int64
    yres::Int64
end
DummyFilm() = DummyFilm(10,10)

@testset "Orthographic camera" begin
    @testset "Orthographic function" begin
        T = orthographic(0.0, 1.0)
        @test T.M ≈ eye(4)
        @test T.MInv ≈ eye(4)
        T = orthographic(1.0, 3.0)
        @test diag(T.M) ≈ [1.0, 1.0, 0.5, 1.0]
        @test T.M[:,4] ≈ [0.0, 0.0, -0.5, 1.0]
    end
    @testset "Creations" begin
        T = Transformation()
        window = [-1.0, 1.0, -1.0, 1.0]
        C = OrthographicCamera(T, window, 0.0, 0.0, DummyFilm())
        sample = CameraSample(5.0, 5.0, 0.5, 0.5)
        R = generate_ray(C, sample)
        @test R.origin ≈ Point3(0.0, 0.0, 0.0)
        @test R.direction ≈ Vector3(0.0, 0.0, 1.0)
        @test R.tmin == 0.0
        @test R.tmax == Inf
        @test R.depth == 1
    end
    @testset "Transformations" begin
        T = Transformation()
        window = [-1.0, 1.0, -1.0, 1.0]
        C = OrthographicCamera(T, window, 0.0, 0.0, DummyFilm())
        @test C.camera_to_world.M ≈ eye(4)
        @test C.raster_to_camera(Point3(0)) ≈ Point3(-1, -1, 0)
        @test C.screen_to_camera(Point3(0)) ≈ Point3(0)
    end
end


@testset "ImageFilm" begin
    @testset "Filter table" begin
        tbl = make_filtertable(BoxFilter(1.0, 1.0))
        @test tbl ≈ ones(16,16)
    end
    @testset "Adding sample" begin
        F = ImageFilm(10, 10, BoxFilter(1.0, 1.0))
        
    end
end

