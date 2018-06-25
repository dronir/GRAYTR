


# world_bounds(P::Primitive) returns BoundingBox
# obj_bounds(P::Primitive) return BoundingBox
# can_intersect(P::Primitive) returns Bool
# intersect(r::Ray, P::Primitive) return Intersection or nothing
# intersectP(r::Ray, P::Primitive)
# refine(P::Primitive)
# BDRF(P::Primitive)

struct DummyShape <: Shape
end

struct GeometricPrimitive{T<:Shape} <: Primitive
    shape::T
    material::Material
    light::Union{AreaLight,Void}
    id::Int64
end

intersectP(R::Ray, P::GeometricPrimitive) = intersectP(R, P.shape)
can_intersect(P::GeometricPrimitive) = can_intersect(P.shape)
world_bounds(P::GeometricPrimitive) = world_bounds(P.shape)
obj_bounds(P::GeometricPrimitive) = obj_bounds(P.shape)

(T::Transformation)(P::GeometricPrimitive) = GeometricPrimitive(T(P.shape), P.material, P.light, P.id)



struct Intersection{P<:Primitive}
    target::P
    geometry::DifferentialGeometry
    obj_to_world::Transformation
    world_to_obj::Transformation
    ray_eps::Float64
    tmin::Float64
    shape_id::Int64
    primitive_id::Int64
end

function intersect(R::Ray, P::GeometricPrimitive)
    dg, tmin, reps = shape_intersect(R, P.shape)
    if dg == nothing
        return nothing
    end
    return Intersection(
        P, dg, P.shape.obj_to_world, P.shape.world_to_obj, reps, tmin, P.shape.id, P.id
    )
end

function get_BSDF(P::GeometricPrimitive, dg::DifferentialGeometry, obj_to_world::Transformation)
    geom = dg #get_shading_geometry(P.shape, dg, obj_to_world)
    return get_BSDF(P.material, dg, geom)
end

get_BSDF(I::Intersection) = get_BSDF(I.target, I.geometry, I.obj_to_world)



# Refine a primitive repeatedly until all primitives are intersectable
function fully_refine!{T<:Primitive}(P::Primitive, out::Array{T,1})
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


