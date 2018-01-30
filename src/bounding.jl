

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
    return BoundingBox(pMin, pMax, contents)
end

# Get the bounding box of a group of points
function BoundingBox(contents::Array{Point3,1})
    pMin = Point3(Inf, Inf, Inf)
    pMax = Point3(-Inf, -Inf, -Inf)
    for p in contents
        pMin = min(pMin, p)
        pMax = max(pMax, p)
    end
    return BoundingBox(pMin, pMax, contents)
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
