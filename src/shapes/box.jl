
# Box type and related methods

struct Box <: Shape
    id::Int64
    obj_to_world::Transformation
    world_to_obj::Transformation
end

Box() = Box(1, Transformation(), Transformation())
Box(T::Transformation) = Box(1, T, inv(T))
Box(id::Int64, T::Transformation) = Box(id, T, inv(T))

can_intersect(b::Box) = true


function (T::Transformation)(B::Box)
    T2 = T * B.obj_to_world
    Box(B.id, T2, inv(T2))
end

obj_bounds(b::Box) = BoundingBox(Point3(-1,-1,-1), Point3(1,1,1))

function world_bounds(B::Box)
    B.obj_to_world(obj_bounds(B))
end



function shape_intersect(r::Ray, B::Box)
    ray = B.world_to_obj(r)

    t0 = ray.tmin
    t1 = ray.tmax
    ibest = -1
    for i = 1:3
        tNear = (-1.0 - ray.origin[i]) / ray.direction[i]
        tFar = (1.0 - ray.origin[i]) / ray.direction[i]
        tNear, tFar = min(tNear, tFar), max(tNear, tFar)
        if tNear > t0 
            t0 = tNear
            ibest = i
        end
        t1 = tFar < t1 ? tFar : t1
        if t0 > t1 
            # Cannot hit the box
            return nothing, NaN, NaN
        end 
    end
    # We hit the box and (t0, t1) are the distances where this happens.
    # These are also already guaranteed to be within (ray.tmin, ray.tmax)
    
    t = t0
    P = ray(t)
    if ibest == 1
        n = Normal3(1,0,0)
        s = Vector3(0,1,0)
    elseif ibest == 2
        n = Normal3(0,1,0)
        s = Vector3(0,0,1)
    elseif ibest == 3
        n = Normal3(0,0,1)
        s = Vector3(1,0,0)
    else
        error("Weirdness in box intersection.")
    end
    
    n = flip_normal(n, ray.direction)
    DG = DifferentialGeometry(B.obj_to_world(P), B.obj_to_world(n), B.obj_to_world(s))
    return DG, t, 5e-4 * t
end



# Quick intersect function
function intersectP(r::Ray, B::Box)
    ray = B.world_to_obj(r)

    t0 = ray.tmin
    t1 = ray.tmax
    for i = 1:3
        tNear = (-1.0 - ray.origin[i]) / ray.direction[i]
        tFar = (1.0 - ray.origin[i]) / ray.direction[i]
        tNear, tFar = min(tNear, tFar), max(tNear, tFar)
        t0 = tNear > t0 ? tNear : t0
        t1 = tFar < t1 ? tFar : t1
        if t0 > t1 
            # Cannot hit the box
            return false
        end 
    end
    return true
end
