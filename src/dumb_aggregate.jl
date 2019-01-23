

"""

    DumbAggregate{P<:Primitive} <: Aggregate

A dumb aggregate is the simplest way to store primitives: it's just a list and all the
primitives are checked for intersection one by one. It's mostly used for testing purposes,
as the `BVHAccelerator` will be faster for any reasonable number of primitives.

"""
struct DumbAggregate{P<:Primitive} <: Aggregate
    primitives::Array{P,1}
    bounds::BoundingBox
end


"""

    DumbAggregate(primitives::Array{P,1})

Construct a `DumbAggregate` from a list of primitives, computing the bounding box.

"""
function DumbAggregate(primitives::Array{P,1}) where P<:Primitive
    BB = BoundingBox([world_bounds(p) for p in primitives])
    return DumbAggregate(primitives, BB)
end


"""
    world_bounds(B::DumbAggregate)

Return the bounding box containing all the primitives in the `DumbAggregate`.

"""
world_bounds(A::DumbAggregate) = A.bounds



"""
    intersect(ray::Ray, A::DumbAggregate)

Find nearest intersection, if any, between `ray` and a `DumbAggregate`. 

Returns `Intersection` or `nothing`.

"""
function intersect(ray::Ray, A::DumbAggregate)
    best_tmin = Inf
    best_isect = nothing
    if !intersectP(ray, A.bounds)
        return nothing
    end
    for primitive in A.primitives
        best_isect = update_isect(best_isect, intersect(ray, primitive))
    end
    return best_isect
end



"""
    intersectP(ray::Ray, A::DumbAggregate)

Returns `true` if the given ray intersects any primitive in the `DumbAggregate`, else
returns false. This is much faster to do than [`intersect`](@ref), since it's enough to
find any intersection, not specifically the nearest one.

Possibly the most often called function in the entire ray-tracer.

"""
function intersectP(R::Ray, A::DumbAggregate)
    if !intersectP(R, A.bounds)
        return false
    end
    for prim = A.primitives
        if intersectP(R, prim)
            return true
        end
    end
    return false
end
