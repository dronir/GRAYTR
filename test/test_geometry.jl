
using LinearAlgebra
using Statistics
using StaticArrays

@testset "Geometry" begin

for type_under_test in [Vector3, Point3, Normal3]
    typename = string(type_under_test)
    @testset "Common tests: $typename" begin
        X = type_under_test(1.0, 0.0, 0.0)
        Y = type_under_test(0.0, 1.0, 0.0)
        Z = type_under_test(0.0, 0.0, 1.0)
        @testset "Creating" begin
            v = type_under_test(1.0, 2.0, 3.0)
            @test v.x ≈ 1.0
            @test v.y ≈ 2.0
            @test v.z ≈ 3.0
            v = type_under_test(1,2,3)
            @test v.x ≈ 1.0
            @test v.y ≈ 2.0
            @test v.z ≈ 3.0
            v = type_under_test(1//2, 3//4, -2//1)
            @test v.x ≈ 0.5
            @test v.y ≈ 0.75
            @test v.z ≈ -2.0
            v = type_under_test(1.0)
            @test v.x ≈ 1.0
            @test v.y ≈ 1.0
            @test v.x ≈ 1.0
            v = type_under_test(Vector3(1, 2, 3))
            @test v.x ≈ 1.0
            @test v.y ≈ 2.0
            @test v.z ≈ 3.0
        end
        
        @testset "Arithmetic" begin
            v = type_under_test(1, 2, 3)
            @test 2.0 * v ≈ type_under_test(2, 4, 6)
            @test v * 2.0 ≈ type_under_test(2, 4, 6)
            @test v/2 ≈ type_under_test(0.5, 1.0, 1.5)
        end
        
        @testset "Normalizing" begin
            v = type_under_test(2, 0, 0)
            @test normalize(v) ≈ type_under_test(1, 0, 0)
            v = type_under_test(0, -2, 0)
            @test normalize(v) ≈ type_under_test(0, -1, 0)
            v = type_under_test(3, -4, 0)
            @test normalize(v) ≈ type_under_test(3//5, -4//5, 0)
        end
        
        @testset "Dot product" begin
            v = type_under_test(0.5, 1.0, -2.0)
            @test dot(v, X) ≈ 0.5
            @test dot(v, Y) ≈ 1.0
            @test dot(v, Z) ≈ -2.0
            
            u = type_under_test(2, 1, 1)
            @test dot(v, u) ≈ 0.0
        end
        
        @testset "Min/max" begin
            v = type_under_test(-1, 0, 2)
            u = type_under_test(0, -1, 0)
            @test min(v,u) ≈ type_under_test(-1, -1, 0)
            @test max(u,v) ≈ type_under_test(0, 0, 2)
        end
        
        @testset "Iteration" begin
            V = type_under_test(1, 2, 3)
            A = [x for x in V]
            @test A ≈ [1.0, 2.0, 3.0]
        end
        
        @testset "Cross product" begin
            X = type_under_test(1,0,0)
            Y = type_under_test(0,1,0)
            Z = type_under_test(0,0,1)
            @test cross(X, Y) ≈ Vector3(0, 0, 1)
            @test cross(X, Z) ≈ Vector3(0, -1, 0)
            @test cross(Y, X) ≈ Vector3(0, 0, -1)
            @test cross(Y, Z) ≈ Vector3(1, 0, 0)
            @test cross(Z, X) ≈ Vector3(0, 1, 0)
            @test cross(Z, Y) ≈ Vector3(-1, 0, 0)
        end
        
        @testset "Various functions" begin
            X = type_under_test(1, 2, 3)
            @test size(X) == 3
            @test length(X) == 3
            @test !isnan(X)
            @test isnan(type_under_test(NaN, 2, 3))
            @test isnan(type_under_test(1, NaN, 3))
            @test isnan(type_under_test(1, 2, NaN))
        end
        
        @testset "Conversions" begin
            T = type_under_test(1, 2, 3)
            @test convert(Array{Float64, 1}, T) ≈ [1.0, 2.0, 3.0]
            
            V = Vector3(1, 2, 3)
            P = Point3(1, 4, 9)
            N = Normal3(1, 8, 27)
            @test convert(type_under_test, V) ≈ type_under_test(1, 2, 3)
            @test convert(type_under_test, P) ≈ type_under_test(1, 4, 9)
            @test convert(type_under_test, N) ≈ type_under_test(1, 8, 27)
        end
    end
end

@testset "Type-specific functions" begin
    P1 = Point3(2, 0, -4)
    P2 = Point3(0, 0, -2)
    @test mean(P1, P2) ≈ Point3(1, 0, -3)
    
    V1 = Vector3(1, 2, 3)
    V2 = Vector3(1, 0, -1)
    @test V1 + V2 ≈ Vector3(2, 2, 2)
    @test V1 - V2 ≈ Vector3(0, 2, 4)
end



@testset "Transformations" begin
    @testset "Identity transformation" begin
        T = Transformation()
        @test T.M ≈ Matrix{Float64}(I, 4, 4)
        @test T.MInv ≈ Matrix{Float64}(I, 4, 4)
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
    @testset "Miscellaneous" begin
        v = Vector3(1/sqrt(3), 1/sqrt(3), 1/sqrt(3))
        for direction in [X_AXIS, Y_AXIS, Z_AXIS, v]
            T = GRAYTR.rotate_z_to(direction)
            @test isapprox(T(Z_AXIS), direction ; atol=1e-16)
        end
    end
end


end # testset Geometry

