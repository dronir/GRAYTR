
using Base.Test
include("geometry.jl")
using Geometry

function test_vector3()
    @testset "Vector3" begin
        X = Vector3(1.0, 0.0, 0.0)
        Y = Vector3(0.0, 1.0, 0.0)
        Z = Vector3(0.0, 0.0, 1.0)
        @testset "Creating" begin
            v = Vector3(1.0, 2.0, 3.0)
            @test v.x ≈ 1
            @test v.y ≈ 2
            @test v.z ≈ 3
            v = Vector3(1,2,3)
            @test v.x ≈ 1
            @test v.y ≈ 2
            @test v.z ≈ 3
            v = Vector3(1//2, 3//4, -2//1)
            @test v.x ≈ 0.5
            @test v.y ≈ 0.75
            @test v.z ≈ -2.0
        end
        
        @testset "Normalizing" begin
            v = Vector3(2, 0, 0)
            @test normalize(v) ≈ Vector3(1, 0, 0)
            v = Vector3(0, -2, 0)
            @test normalize(v) ≈ Vector3(0, -1, 0)
            v = Vector3(3, -4, 0)
            @test normalize(v) ≈ Vector3(3//5, -4//5, 0)
        end
        
        @testset "Dot product" begin
            v = Vector3(0.5, 1.0, -2.0)
            @test dot(v, X) ≈ 0.5
            @test dot(v, Y) ≈ 1.0
            @test dot(v, Z) ≈ -2.0
            
            u = Vector3(2, 1, 1)
            @test dot(v, u) ≈ 0.0
        end
        
        @testset "Min/max" begin
            v = Vector3(-1, 0, 2)
            u = Vector3(0, -1, 0)
            @test min(v,u) ≈ Vector3(-1, -1, 0)
            @test max(u,v) ≈ Vector3(0, 0, 2)
        end
    end
end

test_vector3()
