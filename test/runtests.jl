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

# Tests of materials.jl
include("test_materials.jl")

# Tests of samplers.jl
include("test_samplers.jl")

# Tests of cameras.jl
include("test_cameras.jl")

# Full integration test for image rendering.
#include("test_full_computation.jl")

