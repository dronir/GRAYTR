using GRAYTR
using Test

# Geometry package tests
include("test_geometry.jl")

# Utility function tests
include("test_utility.jl")

# Tests of rays.jl
include("test_rays.jl")

# Tests of bounding.jl
include("test_bounding.jl")

# Tests of shapes
include("test_shapes.jl")

# Tests of spectra
include("test_spectrum.jl")

# Tests of brdf.jl
include("test_brdf.jl")

# Tests of geometric primitives
include("test_primitives.jl")

# Tests of dummy aggregate
include("test_dumbaggregate.jl")

# Tests of samplers.jl
include("test_samplers.jl")

# Tests of filters.jl
include("test_filters.jl")

# Tests of cameras.jl
include("test_cameras.jl")

# Full integration test for image rendering.
include("test_full_computation.jl")

# Full radiation pressure test computation
include("test_pressure.jl")
