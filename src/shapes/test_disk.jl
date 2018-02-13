include("../abstracts.jl")
include("../geometry.jl")
using Geometry
using Base.Test
include("../rays.jl")
include("../bounding.jl")

include("../diffgeom.jl")
include("disk.jl")

@testset "Ray intersection" begin
    C = Disk(Transformation())
    for phi = linspace(0.0, 2pi, 17)
        for sgn = [-1, 1]
            R = Ray(Point3(0.5 * cos(phi), 0.5*sin(phi), sgn*3.0), Vector3(0.0, 0.0, -sgn))
            @test intersectP(R, C)
            maybe_isect, t, reps = intersect(R, C)
            @test !isnull(maybe_isect)
            isect = get(maybe_isect)
            @test t ≈ 3.0
            @test isect.n ≈ Normal3(0.0, 0.0, 1.0)
        end 
    end
end
