

# A bounding box is defined by its extreme corners.
struct BoundingBox
    pMin::Point3
    pMax::Point3
end

# Get the bounding box of a group of bounding boxes
function BoundingBox(boxes::Array{BoundingBox,1})
    pMin = Point3(Inf, Inf, Inf)
    pMax = Point3(-Inf, -Inf, -Inf)
    for bb in boxes
        pMin = min(pMin, bb.pMin)
        pMax = max(pMax, bb.pMax)
    end
    return BoundingBox(pMin, pMax)
end

# Get the bounding box of a group of points
function BoundingBox(contents::Array{Point3,1})
    pMin = Point3(Inf, Inf, Inf)
    pMax = Point3(-Inf, -Inf, -Inf)
    for p in contents
        pMin = min(pMin, p)
        pMax = max(pMax, p)
    end
    return BoundingBox(pMin, pMax)
end

# Apply a transformation on a bounding box (probably not in the most efficient way)
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
            return Nullable{Tuple{Float64,Float64}}()
        end 
    end
    return Nullable((t0, t1))
end

intersectP(R::Ray, BB::BoundingBox) = !isnull(intersect(R, BB))

area(BB::BoundingBox) = 2 * ((BB.pMax[1] - BB.pMin[1]) * (BB.pMax[2] - BB.pMin[2])
                           + (BB.pMax[2] - BB.pMin[2]) * (BB.pMax[3] - BB.pMin[3])
                           + (BB.pMax[3] - BB.pMin[3]) * (BB.pMax[1] - BB.pMin[1]))

