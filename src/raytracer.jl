
module RayTracer


include("geometry.jl")
using .Geometry

abstract type Camera end

struct RayIntersection
    point::Point3
    normal::Normal3
    direction::Vector3
    target::SceneObject
end

## Rays
include("rays.jl")

## Bounding volumes
include("bounding.jl")




## Scene objects
##

struct Box <: SceneObject
    pMin::Point3
    pMax::Point3
end

BoundingBox(B::Box) = BoundingBox(B.pMin, B.pMax, [B])

struct Triangle <: SceneObject
    a::Point3
    b::Point3
    c::Point3
    normal::Normal3
end

function Triangle(a,b,c)
    normal = 0.5 * cross(b-a, c-a)
    return Triangle(a,b,c,normal)
end

struct Sphere <: SceneObject
    center::Point3
    radius::Float64
end



end # module
