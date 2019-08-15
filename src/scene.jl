
export Scene


"""

    Scene

The `Scene` type contains the root node of the object hierarchy, a list of the light
sources in the scene, and the global bounding box.

"""
struct Scene{T<:Aggregate, S<:Spectrum}
    aggregate::T
    lights::Array{LightSource,1}
    bounds::BoundingBox
    background::S
end

Scene(agg::Aggregate, lights::Array{LightSource,1}) = Scene(agg, lights, world_bounds(agg), nolight)
Scene(agg::Aggregate, lights::Array{LightSource,1}, background::Spectrum) = Scene(agg, lights, world_bounds(agg), background)

function Scene(primitives::Array{GeometricPrimitive}, lights::Array{LightSource,1})
    return Scene(BVHAccelerator(primitives), lights)
end

function Scene(primitives::Array{GeometricPrimitive}, lights::Array{LightSource,1}, background::Spectrum)
    return Scene(BVHAccelerator(primitives), lights, background)
end


"""
    world_bounds(S::Scene)

Get the bounding box of the entire scene (this just asks for the bounding box of the
scene's acceleration aggregate.)

"""
world_bounds(S::Scene) = Scene.bounds


"""
    intersect(r::Ray, S::Scene)

Find intersection of a ray in the scene. Just delegates the call to the scene's
acceleration aggregate's `intersect` method. Returns an `Intersection`.

"""
intersect(r::Ray, S::Scene) = intersect(r, S.aggregate)


"""
    intersectP(r::Ray, S::Scene)

A quicker true/false ray intersect function for a scene. This just delegates the call to
the scene's acceleration aggregate's `intersectP`.

"""
intersectP(r::Ray, S::Scene) = intersectP(r, S.aggregate)
