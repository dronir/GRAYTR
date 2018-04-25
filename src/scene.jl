
# The Scene type contains the root node of the object hierarchy, a list of the light sources
# in the scene, and the global bounding box.

struct Scene{T<:Aggregate}
    aggregate::T
    lights::Array{LightSource,1}
    bounds::BoundingBox
end

Scene(agg::Aggregate, lights::Array{LightSource,1}) = Scene(agg, lights, world_bounds(agg))

world_bounds(S::Scene) = Scene.bounds
intersect(r::Ray, S::Scene) = intersect(r, S.aggregate)
intersectP(r::Ray, S::Scene) = intersectP(r, S.aggregate)
