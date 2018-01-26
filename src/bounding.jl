
struct BoundingBox <: BoundingGeometry
    pMin::Point3
    pMax::Point3
    contents::Array{Intersectable,1}
end

struct BoundingSphere <: BoundingGeometry
    center::Point3
    radius::Float64
    contents::Array{Intersectable,1}
end

function BoundingBox(contents::Array{T,1}) where {T<:Intersectable} 
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

function BoundingBox(Sph::BoundingSphere)
    pMin = Sph.center - Sph.radius
    pMax = Sph.center + Sph.radius
    return BoundingBox(pMin, pMax, [Sph])
end

BoundingSphere(Sph::BoundingSphere) = Sph

function BoundingSphere(BB::BoundingBox)
    center = (BB.pMin + BB.pMax) / 2
    radius = norm(BB.pMax - center)
    return BoundingSphere(center, radius, [BB])
end