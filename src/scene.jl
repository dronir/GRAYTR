
struct Scene
    aggregate::Primitive
    lights::Array{LightSource,1}
    bounds::BoundingBox
end

Scene(agg::Primitive, lights::Array{LightSource,1}) = Scene(agg, lights, world_bounds(agg))

world_bounds(S::Scene) = Scene.bounds
intersect(r::Ray, S::Scene) = intersect(r, S.aggregate)
intersectP(r::Ray, S::Scene) = intersectP(r, S.aggregate)
