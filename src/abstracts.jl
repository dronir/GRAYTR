# All the abstract types are declared in this file. Some type definitions depend
# on various abstract types in ways that make keeping each abstract type in the
# most relevant file impossible. 

# E.g. funcions for Scene depends on LightSource and Primitive, Primitive
# depends on AreaLight, and functions for light sources depend on Scene.

# For primitives.jl
abstract type Primitive end

# For bsdr.jl
abstract type BxDF end

# For shape.jl
abstract type Shape end

# For materials.jl
abstract type Material end

# For spectrum.jl
abstract type Spectrum end

# For samplers.jl
abstract type Sampler end
abstract type Sample end

# For renderer.jl
abstract type Renderer end

# For integrator.jl
abstract type Integrator end
abstract type SurfaceIntegrator <: Integrator end
abstract type VolumeIntegrator <: Integrator end

# For lights.jl
abstract type LightSource end

# For cameras.jl
abstract type Camera end
abstract type Film end
abstract type Filter end
