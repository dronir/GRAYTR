

# Interface of a Shape:
#  can_intersect(s::Shape), returns Bool 
#  area(s::Shape), returns Float64
#  obj_bounds(s::Shape), returns BoundingBox
#  world_bounds(s::Shape), returns BoundingBox
#  intersect(r::Ray, s::Shape) returns (Nullable{DifferentialGeometry}, tMin, rayEps)
#  intersectP(r::Ray, s::Shape) returns Bool

import Base.intersect

include("shapes/sphere.jl")




