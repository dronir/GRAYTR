
# This code to be inserted into raytracer.jl with include("rays.jl")

struct Ray
    origin::Point3
    direction::Vector3
    tmin::Float64
    tmax::Float64
    depth::Int64
end

Ray(origin::Point3, direction::Vector3) = Ray(origin, direction, 0.0, Inf, 1)
Ray(parent::Ray, origin::Point3, direction::Vector3) = Ray(origin, direction, 0.0, Inf, parent.depth+1)

(R::Ray)(t::Float64) = R.origin + t*R.direction

function (transform::Transformation)(R::Ray)
    return Ray(transform(R.origin), transform(R.direction))
end
