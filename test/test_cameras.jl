

struct DummyFilm <: GRAYTR.Film 
    resX::Int64
    resY::Int64
end
DummyFilm() = DummyFilm(10,10)


@testset "Orthographic camera" begin
    @testset "Orthographic function" begin
        T = GRAYTR.orthographic(0.0, 1.0)
        @test T.M ≈ Matrix{Float64}(I, 4, 4)
        @test T.MInv ≈ Matrix{Float64}(I, 4, 4)
        T = GRAYTR.orthographic(1.0, 3.0)
        @test diag(T.M) ≈ [1.0, 1.0, 0.5, 1.0]
        @test T.M[:,4] ≈ [0.0, 0.0, -0.5, 1.0]
    end
    @testset "Creations" begin
        T = Transformation()
        window = [-1.0, 1.0, -1.0, 1.0]
        C = OrthographicCamera(T, window, 0.0, 0.0, DummyFilm())
        sample = GRAYTR.CameraSample(5.0, 5.0, 0.5, 0.5)
        weight, R = GRAYTR.generate_ray(C, sample)
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
        @test C.camera_to_world.M ≈ Matrix{Float64}(I, 4, 4)
        @test C.raster_to_camera(Point3(0)) ≈ Point3(-1, -1, 0)
        @test C.screen_to_camera(Point3(0)) ≈ Point3(0)
    end
end


