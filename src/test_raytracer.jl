
using Base.Test
include("raytracer.jl")
import RayTracer.BoundingBox
import RayTracer.Box
using RayTracer.Geometry

function test_boundingbox()
    B1 = Box([-1.0, -1.0, -1.0], [1.0, 1.0, 1.0])
    B2 = Box([1.0, -1.0, -1.0], [3.0, 1.0, 1.0])
    BB1 = BoundingBox(B1)
    BB2 = BoundingBox(B2)
    BB3 = BoundingBox([B1, B2])
    
    @testset begin
        @test BB1.pMin ≈ Point3(-1, -1, -1)
        @test BB1.pMax ≈ Point3( 1,  1,  1)
        @test BB2.pMin ≈ Point3( 1, -1, -1)
        @test BB2.pMax ≈ Point3( 3,  1,  1)
        @test BB3.pMin ≈ BB1.pMin
        @test BB3.pMax ≈ BB2.pMax
        @test BB3.contents[1] === B1
        @test BB3.contents[2] === B2
    end
    
end

test_boundingbox()
