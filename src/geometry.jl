module Geometry

#using LinearAlgebra
#using Statistics
using StaticArrays

import Statistics.mean, Statistics.normalize
import LinearAlgebra.cross, LinearAlgebra.norm, LinearAlgebra.dot, LinearAlgebra.I
import LinearAlgebra.det

import Base.+, Base.-, Base.*, Base./, Base.isnan, Base.inv
import Base.getindex, Base.min, Base.max
import Base.convert, Base.isapprox
import Base.size, Base.length, Base.iterate

export VectorLike, Vector3, Normal3, Point3
export Transformation, translation, scaling, rotation, look_at, swaps_handedness
export rotate_z_to
export X_AXIS, Y_AXIS, Z_AXIS, NEG_X_AXIS, NEG_Y_AXIS, NEG_Z_AXIS, IDENTITY_TRANSFORM


abstract type VectorLike end

# Three kinds of vectors are represented with different types: points in space, vectors between
# points, and normal vectors of surfaces. The reason is that these behave differently under
# transformations, and efficient implementation is easier if they are separated.
# The most important difference is that normals are associated with a surface, and when the
# surface is transformed, the normal must be transformed in a way that maintains its normality
# to the surface. 

struct Vector3 <: VectorLike
    x::Float64
    y::Float64
    z::Float64
end

struct Point3 <: VectorLike
    x::Float64
    y::Float64
    z::Float64
end

struct Normal3 <: VectorLike
    x::Float64
    y::Float64
    z::Float64
end

const X_AXIS = Vector3(1,0,0)
const Y_AXIS = Vector3(0,1,0)
const Z_AXIS = Vector3(0,0,1)
const NEG_X_AXIS = Vector3(-1,0,0)
const NEG_Y_AXIS = Vector3(0,-1,0)
const NEG_Z_AXIS = Vector3(0,0,-1)

Vector3(x::Real) = Vector3(x,x,x)
Point3(x::Real) = Point3(x,x,x)
Normal3(x::Real) = Normal3(x,x,x)

Vector3(v::Vector3) = Vector3(v.x, v.y, v.z)
Normal3(v::Vector3) = Normal3(v.x, v.y, v.z)
Point3(v::Vector3) = Point3(v.x, v.y, v.z)


size(A::VectorLike) = 3
length(A::VectorLike) = 3

mean(A::Point3, B::Point3) = 0.5 * Point3(A.x+B.x, A.y+B.y, A.z+B.z)
+(A::Point3, B::Normal3) = Point3(A.x+B.x, A.y+B.y, A.z+B.z)
+(A::Point3, B::Vector3) = Point3(A.x+B.x, A.y+B.y, A.z+B.z)
+(B::Vector3, A::Point3) = A+B
+(v::Vector3, u::Vector3) = Vector3(v.x+u.x, v.y+u.y, v.z+u.z)
-(v::Vector3, u::Vector3) = Vector3(v.x-u.x, v.y-u.y, v.z-u.z)
-(v::Point3, u::Point3) = Vector3(v.x-u.x, v.y-u.y, v.z-u.z)
-(A::T) where {T<:VectorLike} = T(-A.x, -A.y, -A.z)

dot(A::VectorLike, B::VectorLike) = A.x*B.x + A.y*B.y + A.z*B.z
cross(A::VectorLike, B::VectorLike) = Vector3(A.y*B.z - A.z*B.y, A.z*B.x - A.x*B.z, A.x*B.y - A.y*B.x)
norm(A::VectorLike) = sqrt(A.x^2 + A.y^2 + A.z^2)
normalize(A::VectorLike) = A / norm(A)

min(A::T, B::T) where {T<:VectorLike} = T(min(A.x,B.x), min(A.y,B.y), min(A.z,B.z))
max(A::T, B::T) where {T<:VectorLike} = T(max(A.x,B.x), max(A.y,B.y), max(A.z,B.z))

function isapprox(A::T, B::T; rtol::Real=sqrt(eps(Float64)), atol::Real=0) where {T<:VectorLike}
    return (isapprox(A.x, B.x; rtol=rtol, atol=atol) &&
    isapprox(A.y, B.y; rtol=rtol, atol=atol) &&
    isapprox(A.z, B.z; rtol=rtol, atol=atol))
end

*(a::Number, v::T) where {T<:VectorLike} = T(a*v.x, a*v.y, a*v.z)
*(v::VectorLike, a::Number) = a*v
/(v::T, a::Number) where {T<:VectorLike} = T(v.x/a, v.y/a, v.z/a)

isnan(v::VectorLike) = isnan(v.x) || isnan(v.y) || isnan(v.z)

getindex(v::VectorLike, i::Integer) = i==1 ? v.x : i==2 ? v.y : i==3 ? v.z : throw(BoundsError())

iterate(V::VectorLike) = (V[1], 2)
iterate(V::VectorLike, state::Integer) = state <= 3 ? (getindex(V, state), state+1) : nothing

convert(::Type{Array{Float64,1}}, v::VectorLike) = [v.x, v.y, v.z]
convert(::Type{Vector3}, x::Array{T,1}) where {T<:Number} = Vector3(x...)
convert(::Type{Normal3}, x::Array{T,1}) where {T<:Number} = Normal3(x...)
convert(::Type{Point3}, x::Array{T,1}) where {T<:Number} = Point3(x...)
convert(::Type{Vector3}, v::Normal3) = Vector3(v.x, v.y, v.z)
convert(::Type{Vector3}, v::Point3) = Vector3(v.x, v.y, v.z)
convert(::Type{Normal3}, v::Vector3) = Normal3(v.x, v.y, v.z)
convert(::Type{Normal3}, v::Point3) = Normal3(v.x, v.y, v.z)
convert(::Type{Point3}, v::Vector3) = Point3(v.x, v.y, v.z)
convert(::Type{Point3}, v::Normal3) = Point3(v.x, v.y, v.z)


####################################
# TRANSFORMATIONS


struct Transformation
    M::SArray{Tuple{4,4},Float64,2,16}
    MInv::SArray{Tuple{4,4},Float64,2,16}
end

const EYE = @SMatrix [1.0 0.0 0.0 0.0 ; 0.0 1.0 0.0 0.0 ; 0.0 0.0 1.0 0.0 ; 0.0 0.0 0.0 1.0]
Transformation() = Transformation(EYE, EYE)
Transformation(M::Array{Float64,2}) = Transformation(SMatrix{4,4}(M), SMatrix{4,4}(inv(M)))
Transformation(M::Array{T,2}) where {T<:Real} = Transformation(convert(SArray{4,4}, M))
inv(T::Transformation) = Transformation(T.MInv, T.M)
*(T::Transformation, U::Transformation) = Transformation(T.M*U.M, U.MInv*T.MInv)
# TODO: check maths above

const IDENTITY_TRANSFORM = Transformation()

swaps_handedness(T::Transformation) = det(view(T.M, 1:3, 1:3)) < 0.0

translation(delta::Vector3) = translation(delta.x, delta.y, delta.z)
function translation(x::Real, y::Real, z::Real) 
    M = zeros(Float64, 4, 4)
    M[1,4] = x
    M[2,4] = y
    M[3,4] = z
    M[1,1] = 1.0
    M[2,2] = 1.0
    M[3,3] = 1.0
    M[4,4] = 1.0
    return Transformation(M)
end
function translation(delta::Point3)
    M = Matrix{Float64}(1.0I, 4, 4)
    M[1,4] = delta.x
    M[2,4] = delta.y
    M[3,4] = delta.z
    return Transformation(M)
end

scaling(v::Vector3) = scaling(v.x, v.y, v.z)
function scaling(x::Real, y::Real, z::Real)
    M = zeros(4,4)
    M[1,1] = x
    M[2,2] = y
    M[3,3] = z
    M[4,4] = 1.0
    return Transformation(M)
end

function rotation(axis::Vector3, angle::Number)
    a = normalize(axis)
    s = sin(angle)
    c = cos(angle)
    M = zeros(Float64, 4, 4)
    @inbounds M[1,1] = a.x * a.x + (1.0 - a.x * a.x) * c
    @inbounds M[2,1] = a.x * a.y * (1.0 - c) + a.z * s
    @inbounds M[3,1] = a.x * a.z * (1.0 - c) - a.y * s

    @inbounds M[1,2] = a.y * a.x * (1.0 - c) - a.z * s
    @inbounds M[2,2] = a.y * a.y + (1.0 - a.y * a.y) * c
    @inbounds M[3,2] = a.y * a.z * (1.0 - c) + a.x * s

    @inbounds M[1,3] = a.z * a.x * (1.0 - c) + a.y * s
    @inbounds M[2,3] = a.z * a.y * (1.0 - c) - a.x * s
    @inbounds M[3,3] = a.z * a.z + (1.0 - a.z * a.z) * c

    @inbounds M[4,4] = 1.0
    
    return Transformation(M, M')
end


function (T::Transformation)(v::Vector3)
    
    @inbounds return Vector3(
        T.M[1,1]*v.x + T.M[1,2]*v.y + T.M[1,3]*v.z,
        T.M[2,1]*v.x + T.M[2,2]*v.y + T.M[2,3]*v.z,
        T.M[3,1]*v.x + T.M[3,2]*v.y + T.M[3,3]*v.z
    )
end

function (T::Transformation)(v::Point3)
    @inbounds invw = 1.0 / (T.M[4,1]*v.x + T.M[4,2]*v.y + T.M[4,3]*v.z + T.M[4,4])
    @inbounds p = Point3((T.M[1,1]*v.x + T.M[1,2]*v.y + T.M[1,3]*v.z + T.M[1,4]) * invw,
               (T.M[2,1]*v.x + T.M[2,2]*v.y + T.M[2,3]*v.z + T.M[2,4]) * invw,
               (T.M[3,1]*v.x + T.M[3,2]*v.y + T.M[3,3]*v.z + T.M[3,4]) * invw)
    return p
end

function (T::Transformation)(v::Normal3)
    @inbounds return Normal3(
        T.MInv[1,1]*v.x + T.MInv[2,1]*v.y + T.MInv[3,1]*v.z,
        T.MInv[1,2]*v.x + T.MInv[2,2]*v.y + T.MInv[3,2]*v.z,
        T.MInv[1,3]*v.x + T.MInv[2,3]*v.y + T.MInv[3,3]*v.z
    )
end



"""
    rotate_z_to(direction::Vector3)

Return a transformation that rotates the +Z axis to point in the given direction.

"""
function rotate_z_to(direction::Vector3)
    direction = normalize(direction)
    if direction ≈ Z_AXIS
        return Transformation()
    end
    if direction ≈ -Z_AXIS
        return rotation(X_AXIS, π)
    end
    
    c = direction.z
    s = -sqrt(1.0 - c^2)
    M = zeros(Float64, 4, 4)
    
    # TODO: this could probably be optimized a little. It's needed fairly often for some
    # radiation pressure camera models.
    a = normalize(Vector3(direction.y, -direction.x, 0.0))
    
    M[1,1] = a.x * a.x + (1.0 - a.x * a.x) * c
    M[1,2] = a.x * a.y * (1.0 - c) - a.z * s
    M[1,3] = a.x * a.z * (1.0 - c) + a.y * s
    
    M[2,1] = a.x * a.y * (1.0 - c) + a.z * s
    M[2,2] = a.y * a.y + (1.0 - a.y * a.y) * c
    M[2,3] = a.y * a.z * (1.0 - c) - a.x * s
         
    M[3,1] = a.x * a.z * (1.0 - c) - a.y * s
    M[3,2] = a.y * a.z * (1.0 - c) + a.x * s
    M[3,3] = a.z * a.z + (1.0 - a.z * a.z) * c
    
    M[4,4] = 1.0
    
    return Transformation(M, M')
end

end # module
