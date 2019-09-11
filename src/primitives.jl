

export GeometricPrimitive, apply_material


"""
    GeometricPrimitive{T<:Shape, M<:BxDF}

A geometric primitive, a combination of a shape, a material and possibly a light source.
"""
struct GeometricPrimitive{T<:Shape, M<:BxDF} <: Primitive
    shape::T
    material::M
    light::Union{AreaLight,Nothing}
    id::Int64
end


"""
    apply_material(shapes::Array{T}, mat::Material) where T<:Shape

Return an array of GeometricPrimitives

"""
function apply_material(shapes::Array{T}, mat::BxDF) where T<:Shape
    return GeometricPrimitive[GeometricPrimitive(S, mat, nothing, S.id) for S in shapes]
end


"""
    intersectP(R::Ray, P::GeometricPrimitive)

Quick true/false intersection check between a ray and a geometric primitive.
"""
intersectP(R::Ray, P::GeometricPrimitive) = intersectP(R, P.shape)



"""
    can_intersect(P::GeometricPrimitive)

Check whether a geometric primitive can be intersected by a ray.
"""
can_intersect(P::GeometricPrimitive) = can_intersect(P.shape)



"""
    world_bounds(P::GeometricPrimitive)

Get world bounds of a geometric primitive. These are the world bounds of the associated shape.
"""
world_bounds(P::GeometricPrimitive) = world_bounds(P.shape)



"""
    obj_bounds(P::GeometricPrimitive)

Get object bounds of a geometric primitive. These are the object bounds of the associated shape.
"""
obj_bounds(P::GeometricPrimitive) = obj_bounds(P.shape)



"""
    (T::Transformation)(P::GeometricPrimitive)

Application of a `Transformation` on a `GeometricPrimitive` passes the transformation
on to the shape and returns a new `GeometricPrimitive` with the transformed shape.
"""
(T::Transformation)(P::GeometricPrimitive) = GeometricPrimitive(T(P.shape), P.material, P.light, P.id)



"""
    Intersection{M<:BxDF}

A structure storing information about a ray-scene intersection
"""
struct Intersection{M<:BxDF}
    material::M
    geometry::DifferentialGeometry
    tmin::Float64
end


"""
    update_isect(isect::Intersection, best_isect::Intersection)

Compare `isect` and `best_isect` and return the one with a smaller `tmin`.

"""
function update_isect(isect::Intersection, best_isect::Intersection)
    if isect.tmin < best_isect.tmin
        return isect
    else
        return best_isect
    end
end

"""
    update_isect(isect::Nothing, best_isect::Intersection)

Always return `best_isect` when `isect` is `nothing`.
"""
update_isect(isect::Nothing, best_isect::Intersection) = best_isect


"""
    update_isect(isect, best_isect::Nothing)

Always return `isect` when `best_isect` is `nothing`.
"""
update_isect(isect::Intersection, best_isect::Nothing) = isect

update_isect(isect::Nothing, best_isect::Nothing) = nothing




"""
    intersect(R::Ray, P::GeometricPrimitive)

Find and get intersection of a ray and a geometric primitive.

Returns an `Intersection` or `nothing`.
"""
function intersect(R::Ray, P::GeometricPrimitive)
    dg, tmin, reps = shape_intersect(R, P.shape)
    if dg == nothing
        return nothing
    end
    return Intersection(
        P.material, dg, tmin
    )
end



"""
    fully_refine!(P::Primitive, out::Array{T,1}) where T <: Primitive

Refine a composite primitive (i.e. split it into its components) and put the
components into the provided array. Repeat until all primitives in the list are
intersectable.

TODO: This should be rewritten to just take a list of primitives and refine all of them.
"""
function fully_refine!(P::Primitive, out::Array{T,1}) where T <: Primitive
    if can_intersect(P)
        push!(out, P)
        return out
    end
    todo = Primitive[P]
    while length(todo) > 0
        piece = pop!(todo)
        if can_intersect(piece)
            push!(out, piece)
        else
            prepend!(todo, refine(piece))
        end
    end
    return out
end


