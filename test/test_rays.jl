
import GRAYTR.Ray


@testset "Rays" begin
    origin = Point3(0,0,0)
    direction = Vector3(1.0, 0.0, 0.0)
    @testset "Creation" begin
        R = Ray(origin, direction)
        @test R.origin === origin
        @test R.direction === direction
        @test R.tmin ≈ 0
        @test R.tmax == Inf
        @test R.depth == 1
    end
    @testset "Ray extension" begin
        R = Ray(origin, direction)
        R(1.0) ≈ Point3(1, 0, 0)
        R(-1.0) ≈ Point3(-1, 0, 0)
        R = Ray(origin, Vector3(2.0, -1.0, 0.0))
        R(2.0) ≈ Point3(4.0, -2.0, 0.0)
    end
end

