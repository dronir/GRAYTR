


################################
# Visibility Tester

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


################################
# Point light

struct PointLight{S<:Spectrum} <: LightSource
    position::Point3
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
    nsamples::Int64
end

# Minimal constructor
function PointLight(L::Spectrum, l2w::Transformation, nsamples::Int64)
    PointLight(l2w(Point3(0)), L, l2w, inv(l2w), nsamples)
end

isdelta(L::PointLight) = true
power(L::PointLight) = 4π * L.intensity

function sample_L(light::PointLight, p::Point3)
    distance = norm(light.position - p)
    L = light.intensity / distance^2
    wi = (light.position - p) / distance
    pdf = 1.0
    vis = VisibilityTester(p, light.position, 2e-5, 0.0)
    return L, wi, pdf, vis
end

direct(L::PointLight) = true
background(L::PointLight) = NoLight()


################################
# Distance light source

struct DistantLight{S<:Spectrum} <: LightSource
    direction::Vector3
    intensity::S
    light_to_world::Transformation
    world_to_light::Transformation
    nsamples::Int64
end

function DistantLight(dir::Vector3, L::Spectrum, l2w::Transformation, nsamples::Int64)
    DistantLight(normalize(l2w(dir)), L, l2w, inv(l2w), nsamples)
end

function sample_L(light::DistantLight, p::Point3)
    return light.intensity, light.direction, 1.0, VisibilityTester(p, light.direction, 2e-5)
end

direct(L::DistantLight) = true
background(L::DistantLight) = NoLight()


################################
# Distance light source

struct Background{S<:Spectrum} <: LightSource
    intensity::S
end

sample_L(light::Background, p::Point3) = NoLight()

direct(L::Background) = false
background(L::Background) = L.intensity



################################
# Area Light (TODO)

struct AreaLight <: LightSource end

