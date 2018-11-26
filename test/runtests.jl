using GRAYTR
using Test

# Geometry package tests
include("test_geometry.jl")

# Tests of rays.jl
include("test_rays.jl")

# Tests of materials.jl
include("test_materials.jl")

# Tests of samplers.jl
include("test_samplers.jl")

# Full integration test for image rendering.
#include("test_full_computation.jl")

