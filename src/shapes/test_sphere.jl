include("../geometry.jl")
using Geometry
using Base.Test
include("../bounding.jl")
include("../rays.jl")

abstract type Shape end

include("../diffgeom.jl")
include("sphere.jl")



function test_basic()
    T = Transformation()
    @testset "Creating" begin
        Sph = Sphere(1, 2.0, false, T, inv(T))
        @test Sph.id == 1
        @test Sph.radius ≈ 2.0
        @test !Sph.inverted
        @test Sph.obj_to_world.M ≈ eye(4)
        @test Sph.world_to_obj.M ≈ eye(4)
    end
    @testset "Short constructor 1" begin
        Sph = Sphere(1, 3.0, T)
        @test Sph.id == 1
        @test Sph.radius ≈ 3.0
        @test !Sph.inverted
        @test Sph.obj_to_world.M ≈ eye(4)
        @test Sph.world_to_obj.M ≈ eye(4)
    end
    @testset "Short constructor 2" begin
        Sph = Sphere(1, 3.0, true, T)
        @test Sph.id == 1
        @test Sph.radius ≈ 3.0
        @test Sph.inverted
        @test Sph.obj_to_world.M ≈ eye(4)
        @test Sph.world_to_obj.M ≈ eye(4)
    end
    @testset "Simple functions" begin
        Sph = Sphere(1, 2.0, T)
        @test area(Sph) ≈ 16pi
        @test can_intersect(Sph)
    end
    @testset "Bounding box" begin
        @testset "Regular sphere" begin
            Sph = Sphere(1, 1.0, T)
            BB = object_bounds(Sph)
            @test BB.pMin ≈ Point3(-1, -1, -1)
            @test BB.pMax ≈ Point3(1, 1, 1)
            BB = world_bounds(Sph)
            @test BB.pMin ≈ Point3(-1, -1, -1)
            @test BB.pMax ≈ Point3(1, 1, 1)
        end
        @testset "Flattened sphere" begin
            flatten = scaling(Vector3(1.0, 1.0, 0.5))
            Sph = Sphere(1, 1.0, flatten)
            BB = world_bounds(Sph)
            @test BB.pMin ≈ Point3(-1, -1, -0.5)
            @test BB.pMax ≈ Point3(1, 1, 0.5)
        end
        @testset "Shifted sphere" begin
            flatten = translation(Vector3(1.0, 0.0, 0.0))
            Sph = Sphere(1, 1.0, flatten)
            BB = world_bounds(Sph)
            @test BB.pMin ≈ Point3(0.0, -1.0, -1.0)
            @test BB.pMax ≈ Point3(2.0, 1.0, 1.0)
        end    
    end
end

function test_intersections()
    @testset "Simple true/false ray intersection" begin
        @testset "Regular sphere" begin
            Sph = Sphere(1, 1.0, Transformation())
            ray = Ray(Point3(2, 0, 0), Vector3(-1, 0, 0))
            @test intersectP(ray, Sph)
            ray = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
            @test !intersectP(ray, Sph)
        end
        @testset "Shifted sphere" begin
            shift = translation(1, 0, 0)
            Sph = Sphere(1, 1.0, shift)
            ray = Ray(Point3(10, 0, 0), Vector3(-1, 0, 0))
            @test intersectP(ray, Sph)
            ray = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
            @test !intersectP(ray, Sph)
        end
        @testset "Stretched sphere" begin
            shift = scaling(2, 1, 1)
            Sph = Sphere(1, 1.0, shift)
            ray = Ray(Point3(10, 0, 0), Vector3(-1, 0, 0))
            @test intersectP(ray, Sph)
            ray = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
            @test !intersectP(ray, Sph)
        end
    end
    @testset "Complicated ray intersection" begin
        @testset "Regular sphere" begin
            Sph = Sphere(1, 1.0, Transformation())
            ray = Ray(Point3(2, 0, 0), Vector3(-1, 0, 0))
            DG_maybe, t, ray_eps = intersect(ray, Sph)
            @test !isnull(DG_maybe)
            @test t ≈ 1.0
            
            ray = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
            DG_maybe, t, ray_eps = intersect(ray, Sph)
            @test isnull(DG_maybe)
        end
        @testset "Shifted sphere" begin
            shift = translation(1, 0, 0)
            Sph = Sphere(1, 1.0, shift)
            ray = Ray(Point3(10, 0, 0), Vector3(-1, 0, 0))
            DG_maybe, t, ray_eps = intersect(ray, Sph)
            @test !isnull(DG_maybe)
            @test t ≈ 8.0
            
            ray = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
            DG_maybe, t, ray_eps = intersect(ray, Sph)
            @test isnull(DG_maybe)
        end
        @testset "Stretched sphere" begin
            shift = scaling(2, 1, 1)
            Sph = Sphere(1, 1.0, shift)
            ray = Ray(Point3(10, 0, 0), Vector3(-1, 0, 0))
            DG_maybe, t, ray_eps = intersect(ray, Sph)
            @test !isnull(DG_maybe)
            @test t ≈ 8.0
            
            ray = Ray(Point3(2, 2, 0), Vector3(-1, 0, 0))
            DG_maybe, t, ray_eps = intersect(ray, Sph)
            @test isnull(DG_maybe)
        end
    end
end

test_basic()
test_intersections()