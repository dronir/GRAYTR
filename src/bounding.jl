
abstract type BoundingGeometry end

struct BoundingBox <: BoundingGeometry
    pMin::Point3
    pMax::Point3
    contents::Array{Intersectable,1}
end

function BoundingBox(contents::Array{T,1}) where {T<BoundingGeometry} 
    pMin = Point3(Inf, Inf, Inf)
    pMax = Point3(-Inf, -Inf, -Inf)
    for obj in contents
        bb = BoundingBox(obj)
        pMin = min(pMin, bb.pMin)
        pMax = max(pMax, bb.pMax)
    end
    return BoundingBox(pMin, pMax, contents)
end

BoundingBox(BB::BoundingBox) = BB


