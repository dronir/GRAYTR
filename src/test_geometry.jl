
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
        
        @testset "Iteration" begin
            V = Vector3(1, 2, 3)
            A = [x for x in V]
            @test A ≈ [1.0, 2.0, 3.0]
        end
        
        @testset "Cross product" begin
            X = Vector3(1,0,0)
            Y = Vector3(0,1,0)
            Z = Vector3(0,0,1)
            @test cross(X,Y) ≈ Z
            @test cross(X,Z) ≈ -Y
            @test cross(Y,X) ≈ -Z
            @test cross(Y,Z) ≈ X
            @test cross(Z,X) ≈ Y
            @test cross(Z,Y) ≈ -X
        end
    end
end

function test_transformations()
    @testset "Transformation" begin
        @testset "Basic" begin
            T = Transformation()
            @test T.M ≈ eye(4)
            @test T.MInv ≈ eye(4)
            @test !swaps_handedness(T)
        end
        @testset "Translation" begin
            point = Point3(0)
            T = translation(1, 0, 0)
            @test T(point) ≈ Point3(1.0, 0.0, 0.0)
            U = inv(T)
            @test U(T(point)) ≈ point

            point = Point3(-1.0, 0.0, 2.0)
            T = translation(1, 2, 5)
            @test T(point) ≈ Point3(0.0, 2.0, 7.0)
            U = inv(T)
            @test U(T(point)) ≈ point
        end
        @testset "Rotation" begin
            X = Vector3(1, 0, 0)
            Y = Vector3(0, 1, 0)
            Z = Vector3(0, 0, 1)
            
            point = Point3(1,0,0)
            T = rotation(X, pi/2)
            @test T(point) ≈ Point3(1.0, 0.0, 0.0)
            U = inv(T)
            @test U(T(point)) ≈ point
            
            point = Point3(1.0, 0.0, 0.0)
            T = rotation(Z, pi/2)
            @test isapprox(T(point), Point3(0.0, 1.0, 0.0) ; atol=1e-16)
            U = inv(T)
            @test U(T(point)) ≈ point
            T = rotation(Y, pi/2)
            @test isapprox(T(point), Point3(0.0, 0.0, -1.0) ; atol=1e-16)
            U = inv(T)
            @test U(T(point)) ≈ point
        end
        @testset "Scaling" begin
            P = Point3(-1, 1, 2)
            T = scaling(2, 1, 1)
            @test T(P) ≈ Point3(-2, 1, 2)
            U = inv(T)
            @test U(T(P)) ≈ P
            T = scaling(1, 2, 1)
            @test T(P) ≈ Point3(-1, 2, 2)
            U = inv(T)
            @test U(T(P)) ≈ P
            T = scaling(1, 1, 2)
            @test T(P) ≈ Point3(-1, 1, 4)
            U = inv(T)
            @test U(T(P)) ≈ P
            
        end
    end
    true
end

test_vector3()
test_transformations()
