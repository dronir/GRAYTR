
"""
    BoundingBox(pMin::Point3, pMax::Point3)

A bounding box is defined by its extreme corners `pMin` and `pMax`.

"""
struct BoundingBox <: BoundingVolume
    pMin::Point3
    pMax::Point3
    function BoundingBox(p1::Point3, p2::Point3)
        pMin = min(p1, p2)
        pMax = max(p1, p2)
        return new(pMin, pMax)
    end
end


"""
    BoundingBox(boxes::Array{BoundingBox,1})

Construct the bounding box of a list of bounding boxes.

"""
function BoundingBox(boxes::Array{BoundingBox,1})
    pMin = Point3(Inf, Inf, Inf)
    pMax = Point3(-Inf, -Inf, -Inf)
    for bb in boxes
        pMin = min(pMin, bb.pMin)
        pMax = max(pMax, bb.pMax)
    end
    return BoundingBox(pMin, pMax)
end


"""
    BoundingBox(p1::Point3, p2::Point3, p3::Point3)
    
Construct the bounding box of three points.

"""
function BoundingBox(p1::Point3, p2::Point3, p3::Point3)
    pMin = min(p1, p2, p3)
    pMax = max(p1, p2, p3)
    return BoundingBox(pMin, pMax)
end


"""
    BoundingBox(contents::Array{Point3,1})

Construct the bounding box of a list of points.

"""
function BoundingBox(contents::Array{Point3,1})
    pMin = Point3(Inf, Inf, Inf)
    pMax = Point3(-Inf, -Inf, -Inf)
    for p in contents
        pMin = min(pMin, p)
        pMax = max(pMax, p)
    end
    return BoundingBox(pMin, pMax)
end


"""
    BoundingBox(BB::BoundingBox, p::Point3)

Return a bounding box that covers both the box `BB` and the point `p`.

"""
BoundingBox(BB::BoundingBox, p::Point3) = BoundingBox(min(BB.pMin, p), max(BB.pMax, p))


"""
    BoundingBox(p::Point3, BB::BoundingBox)

Return a bounding box that covers both the box `BB` and the point `p`.

"""
BoundingBox(p::Point3, BB::BoundingBox) = BoundingBox(BB, p)


"""
    BoundingBox(BB1::BoundingBox, BB2::BoundingBox)

Return a bounding box that covers the two bounding boxes given.

"""
BoundingBox(BB1::BoundingBox, BB2::BoundingBox) = BoundingBox(min(BB1.pMin, BB2.pMin), max(BB1.pMax, BB2.pMax))


"""
    (transform::Transformation)(BB::BoundingBox)

Apply a transformation on a bounding box (probably not in the most efficient way).

"""
function (transform::Transformation)(BB::BoundingBox)
    corners = Point3[
        Point3(BB.pMin.x, BB.pMin.y, BB.pMin.z),
        Point3(BB.pMin.x, BB.pMin.y, BB.pMax.z),
        Point3(BB.pMin.x, BB.pMax.y, BB.pMin.z),
        Point3(BB.pMin.x, BB.pMax.y, BB.pMax.z),
        Point3(BB.pMax.x, BB.pMin.y, BB.pMin.z),
        Point3(BB.pMax.x, BB.pMin.y, BB.pMax.z),
        Point3(BB.pMax.x, BB.pMax.y, BB.pMin.z),
        Point3(BB.pMax.x, BB.pMax.y, BB.pMax.z)
    ]
    return BoundingBox([transform(p) for p in corners])
end


"""
    intersect(R::Ray, BB::BoundingBox)

Ray intersection with a bounding box. 

Returns `(hit::Bool, tmin::Float64, tmax::Float64)`, where tmin and tmax are the near and far
hit locations on the box. They will be NaN if hit is false.

"""
function intersect(R::Ray, BB::BoundingBox)
    t0 = R.tmin
    t1 = R.tmax
    for i = 1:3
        tNear = (BB.pMin[i] - R.origin[i]) / R.direction[i]
        tFar = (BB.pMax[i] - R.origin[i]) / R.direction[i]
        tNear, tFar = min(tNear, tFar), max(tNear, tFar)
        t0 = tNear > t0 ? tNear : t0
        t1 = tFar < t1 ? tFar : t1
        if t0 > t1 
            return false, NaN, NaN
        end 
    end
    return true, t0, t1
end


"""
    intersectP(R::Ray, BB::BoundingBox)

True/false intersection test with bounding box.

"""
function intersectP(R::Ray, BB::BoundingBox)
    t0 = R.tmin
    t1 = R.tmax
    for i = 1:3
        invd = 1.0 / R.direction[i]
        tNear = (BB.pMin[i] - R.origin[i]) * invd
        tFar = (BB.pMax[i] - R.origin[i]) * invd
        tNear, tFar = min(tNear, tFar), max(tNear, tFar)
        t0 = tNear > t0 ? tNear : t0
        t1 = tFar < t1 ? tFar : t1
        if t0 > t1 
            return false
        end 
    end
    return true
end

"""
    area(BB::BoundingBox)

Return the surface area of the bounding box.

"""
area(BB::BoundingBox) = 2 * ((BB.pMax[1] - BB.pMin[1]) * (BB.pMax[2] - BB.pMin[2])
                           + (BB.pMax[2] - BB.pMin[2]) * (BB.pMax[3] - BB.pMin[3])
                           + (BB.pMax[3] - BB.pMin[3]) * (BB.pMax[1] - BB.pMin[1]))


"""
    max_extent(BB::BoundingBox)

Get the index of the axis where the bounding box is widest. X is 1, Y is 2 and Z is 3.

"""
function max_extent(BB::BoundingBox)
    dx = BB.pMax.x - BB.pMin.x
    dy = BB.pMax.y - BB.pMin.y
    dz = BB.pMax.z - BB.pMin.z
    if dx > dy && dx > dz
        return 1
    elseif dy > dx && dy > dz
        return 2
    else
        return 3
    end 
end


##### BOUNDING SPHERE


"""
    BoundingSphere(center::Point3, radius::Float64)

A spherical bounding volume.

"""
struct BoundingSphere <: BoundingVolume
    center::Point3
    radius::Float64
end


"""
    BoundingSphere(BB:BoundingBox ; r_eps=1e-15)

Construct a `BoundingSphere` which circumscribes the given `BoundingBox`.

The small value `r_eps` is added to the computed radius to avoid clipping with the corners
of the box.

"""
function BoundingSphere(BB::BoundingBox ; r_eps=1e-15)
    c = mean(BB.pMin, BB.pMax)
    r = norm(c - BB.pMax) + r_eps
    return BoundingSphere(c, r)
end


"""
    intersectP(R::Ray, BS::BoundingSphere)

True/false intersection test between ray and given `BoundingSphere`.

"""
function intersectP(R::Ray, BS::BoundingSphere)
    A = dot(R.direction, R.direction)
    B = 2.0 * dot(R.direction, R.origin)
    C = dot(R.origin, R.origin) - BS.radius^2
    hit, tnear, tfar = quadratic(A, B, C)
    
    if !hit || tnear > R.tmax || tfar < R.tmin || (tnear < R.tmin && tfar > R.tmax)
        return false
    else
        return true
    end
end


@inline cross_section(B::BoundingSphere) = π * B.radius^2
