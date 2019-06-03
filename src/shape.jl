

# Interface of a Shape:
#  can_intersect(s::Shape), returns Bool 
#  area(s::Shape), returns Float64
#  obj_bounds(s::Shape), returns BoundingBox
#  world_bounds(s::Shape), returns BoundingBox
#  intersect(r::Ray, s::Shape) returns (Union{DifferentialGeometry,Void}, tMin, rayEps)
#  intersectP(r::Ray, s::Shape) returns Bool

import Base.intersect

world_to_obj(S::Shape)::Transformation = S.world_to_obj
obj_to_world(S::Shape)::Transformation = S.obj_to_world


include("shapes/sphere.jl")
include("shapes/cylinder.jl")
include("shapes/disk.jl")
include("shapes/triangle.jl")
include("shapes/paraboloid.jl")

#include("shapes/plane.jl")
#include("shapes/box.jl")



