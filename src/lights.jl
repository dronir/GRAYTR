
abstract type LightSource end

struct VisibilityTester
    ray::Ray
end

# VisibilityTester between two points (with epsilons)
function VisibilityTester(p1::Point3, p2::Point3, eps1::Float64, eps2::Float64)
    distance = norm(p2-p1)
    VisibilityTester(Ray(p1, (p2-p1)/distance, eps1, distance*(1.0 - eps2), 1))
end

# VisibilityTester from a point to infinitely far away
function VisibilityTester(origin::Point3, direction::Vector3, eps::Float64)
    VisibilityTester(Ray(origin, direction, eps, Inf, 1))
end

unoccluded(V::VisibilityTester, scene::Scene) = !intersectP(V.ray, scene)


struct PointLight
    position::Point3
    intensity::Spectrum
    light_to_world::Transformation
    world_to_light::Transformation
    nsamples::Int64
end

isdelta(L::PointLight) = true
power(L::PointLight) = 4π * L.intensity

function sample_L(light::PointLight, p::Point3)
    distance = norm(light.position - p)
    L = light.intensity / distance^2
    wi = (light.position - p) / distance
    pdf = 1.0
    vis = VisibilityTester(p, light.position, pEps, 0.0)
    return L, wi, pdf, vis
end

