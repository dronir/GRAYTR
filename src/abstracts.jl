abstract type Primitive end

abstract type BxDF end

abstract type Shape end

abstract type Material end

# For spectrum.jl
abstract type Spectrum end

# For samplers.jl
abstract type Sampler end
abstract type Sample end

# For renderer.jl
abstract type Renderer end
abstract type Task end

# For integrator.jl
abstract type Integrator end
abstract type SurfaceIntegrator <: Integrator end
abstract type VolumeIntegrator <: Integrator end

abstract type LightSource end

