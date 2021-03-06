

@testset "Bounding volumes" begin

@testset "Create bounding box" begin

    p1 = Point3(0, 0, 0)
    p2 = Point3(1, 2, 3)
    p3 = Point3(-3, -2, -1)

    @testset "Degenerate bounding box" begin
        BB = GRAYTR.BoundingBox(p1, p1)
        @test BB.pMin == p1
        @test BB.pMax == p1
        
        BB = GRAYTR.BoundingBox(p2, p2)
        @test BB.pMin == p2
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p3, p3)
        @test BB.pMin == p3
        @test BB.pMax == p3
    end

    @testset "Bounding box from two points" begin
        BB = GRAYTR.BoundingBox(p1, p2)
        @test BB.pMin == p1
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p2, p1)
        @test BB.pMin == p1
        @test BB.pMax == p2
    end
    
    @testset "Bounding box from three points" begin
        BB = GRAYTR.BoundingBox(p1, p2, p3)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p1, p3, p2)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p2, p1, p3)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p2, p3, p1)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p3, p1, p2)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p3, p2, p1)
        @test BB.pMin == p3
        @test BB.pMax == p2
    end
    
    @testset "Bounding box from array of points" begin
        BB = GRAYTR.BoundingBox([p1, p2, p3])
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox([p1, p3, p2])
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox([p2, p1, p3])
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox([p2, p3, p1])
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox([p3, p1, p2])
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox([p3, p2, p1])
        @test BB.pMin == p3
        @test BB.pMax == p2
    end
    
    @testset "Bounding box from two bounding boxes" begin
        BB1 = GRAYTR.BoundingBox(p1, p2)
        BB2 = GRAYTR.BoundingBox(p3, p2)
        
        BB = GRAYTR.BoundingBox(BB1, BB2)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(BB2, BB1)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox([BB1, BB2])
        @test BB.pMin == p3
        @test BB.pMax == p2
    end
    
    @testset "Bounding box from a bounding box and a point" begin
        BB1 = GRAYTR.BoundingBox(p1, p2)
        
        BB = GRAYTR.BoundingBox(BB1, p3)
        @test BB.pMin == p3
        @test BB.pMax == p2
        
        BB = GRAYTR.BoundingBox(p3, BB1)
        @test BB.pMin == p3
        @test BB.pMax == p2
    end
end


@testset "Area" begin
    p1 = Point3(0, 0, 0)
    p2 = Point3(1, 1, 1)
    p3 = Point3(-3, -2, -1)

    BB = GRAYTR.BoundingBox(p1, p2)
    @test GRAYTR.area(BB) ≈ 6.0
    
    BB = GRAYTR.BoundingBox(p1, p1)
    @test GRAYTR.area(BB) ≈ 0.0
    
    BB = GRAYTR.BoundingBox(p1, p3)
    @test GRAYTR.area(BB) ≈ 2*2.0 + 2*3.0 + 2*6.0
end

@testset "Max extent" begin
    p1 = Point3(0, 0, 0)
    p2 = Point3(1, 2, 3)
    p3 = Point3(-3, -2, -1)
    p4 = Point3(1, 3, 2)
    
    BB = GRAYTR.BoundingBox(p1, p2)
    @test GRAYTR.max_extent(BB) == 3
    
    BB = GRAYTR.BoundingBox(p1, p3)
    @test GRAYTR.max_extent(BB) == 1
    
    BB = GRAYTR.BoundingBox(p2, p3)
    @test GRAYTR.max_extent(BB) == 3
    
    BB = GRAYTR.BoundingBox(p1, p1)
    @test GRAYTR.max_extent(BB) == 3
    
    BB = GRAYTR.BoundingBox(p1, p4)
    @test GRAYTR.max_extent(BB) == 2
end


@testset "True/false ray intersection test" begin
    p1 = Point3(0, 0, 0)
    p2 = Point3(2, 2, 1)
    BB = GRAYTR.BoundingBox(p1, p2)
    
    @testset "Ray from inside box" begin
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(1, 0, 0))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, 1, 0))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, 0,  1))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(-1, 0, 0))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, -1, 0))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, 0, -1))
        @test GRAYTR.intersectP(R, BB)
    end
    
    @testset "Ray travelling in -Z direction" begin
        R = GRAYTR.Ray(Point3(0.0, 0.0, 2.0), Vector3(0, 0, -1))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 2.0), Vector3(0, 0, -1))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(2.0, 2.0, 2.0), Vector3(0, 0, -1))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(-1.0, 1.0, 2.0), Vector3(0, 0, -1))
        @test !GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(3.0, 1.0, 2.0), Vector3(0, 0, -1))
        @test !GRAYTR.intersectP(R, BB)
    end
    
    @testset "Ray travelling +X direction" begin
        R = GRAYTR.Ray(Point3(-1.0, 1.0, 0.5), Vector3(1, 0, 0))
        @test GRAYTR.intersectP(R, BB)
        
        R = GRAYTR.Ray(Point3(-1.0, 3.0, 0.5), Vector3(1, 0, 0))
        @test !GRAYTR.intersectP(R, BB)
    end
    
    
end


@testset "Full ray intersection test" begin
    p1 = Point3(0, 0, 0)
    p2 = Point3(2, 2, 1)
    BB = GRAYTR.BoundingBox(p1, p2)
    
    @testset "Ray from inside box" begin
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(1, 0, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 0.0
        @test t1 ≈ 1.0
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, 1, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 0.0
        @test t1 ≈ 1.0
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, 0,  1))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 0.0
        @test t1 ≈ 0.5
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(-1, 0, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 0.0
        @test t1 ≈ 1.0
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, -1, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 0.0
        @test t1 ≈ 1.0
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 0.5), Vector3(0, 0, -1))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 0.0
        @test t1 ≈ 0.5
    end
    
    @testset "Ray travelling in -Z direction" begin
        R = GRAYTR.Ray(Point3(0.0, 0.0, 2.0), Vector3(0, 0, -1))
        isect, t0, t1 =  GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 1.0
        @test t1 ≈ 2.0
        
        R = GRAYTR.Ray(Point3(1.0, 1.0, 2.0), Vector3(0, 0, -1))
        isect, t0, t1 =  GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 1.0
        @test t1 ≈ 2.0
        
        
        R = GRAYTR.Ray(Point3(-1.0, 1.0, 2.0), Vector3(0, 0, -1))
        isect, t0, t1 =  GRAYTR.intersect(R, BB)
        @test !isect
        @test isnan(t0)
        @test isnan(t1)
        
        R = GRAYTR.Ray(Point3(3.0, 1.0, 2.0), Vector3(0, 0, -1))
        isect, t0, t1 =  GRAYTR.intersect(R, BB)
        @test !isect
        @test isnan(t0)
        @test isnan(t1)
    end
    
    @testset "Ray travelling +X direction" begin
        R = GRAYTR.Ray(Point3(-1.0, 1.0, 0.5), Vector3(1, 0, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 1.0
        @test t1 ≈ 3.0
        
        R = GRAYTR.Ray(Point3(-1.0, 3.0, 0.5), Vector3(1, 0, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test !isect
        @test isnan(t0)
        @test isnan(t1)
    end
    
    @testset "Ray travelling -X direction" begin
        R = GRAYTR.Ray(Point3(3.0, 1.0, 0.5), Vector3(-1, 0, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test isect
        @test t0 ≈ 1.0
        @test t1 ≈ 3.0
        
        R = GRAYTR.Ray(Point3(3.0, 3.0, 0.5), Vector3(-1, 0, 0))
        isect, t0, t1 = GRAYTR.intersect(R, BB)
        @test !isect
        @test isnan(t0)
        @test isnan(t1)
    end
    
    
end


@testset "Bounding Sphere" begin
    @testset "Creating" begin
        BSph = GRAYTR.BoundingSphere(Point3(0), 1.0)
        @test BSph.center ≈ Point3(0)
        @test BSph.radius ≈ 1.0
    
        p1 = Point3(0, 0, 0)
        p2 = Point3(1, 2, 4)
    
        d = sqrt(1 + 4 + 16)
    
        BB = GRAYTR.BoundingBox(p1, p2)
    
        BSph = GRAYTR.BoundingSphere(BB)
        @test BSph.center ≈ Point3(0.5, 1.0, 2.0)
        @test isapprox(BSph.radius, d/2 ; atol=1e-14)
    end
    
    @testset "Ray intersection" begin
        BSph = GRAYTR.BoundingSphere(Point3(0), 1.0)
        
        R = GRAYTR.Ray(Point3(0.0, 0.0, 3.0), Vector3(0.0, 0.0, -1.0))
        @test GRAYTR.intersectP(R, BSph)
        
        R = GRAYTR.Ray(Point3(0.0, 0.0, -3.0), Vector3(0.0, 0.0, 1.0))
        @test GRAYTR.intersectP(R, BSph)
        
        R = GRAYTR.Ray(Point3(0.0, 0.0, -3.0), Vector3(0.0, 0.0, -1.0))
        @test !GRAYTR.intersectP(R, BSph)
        
        R = GRAYTR.Ray(Point3(3.0, 0.0, 3.0), Vector3(0.0, 0.0, -1.0))
        @test !GRAYTR.intersectP(R, BSph)
        
    end
    
end


end # main testset
