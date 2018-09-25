module Geometry

using LinearAlgebra
using Statistics

import Statistics.mean, Statistics.normalize
import LinearAlgebra.cross, LinearAlgebra.norm, LinearAlgebra.dot

#Base.dot, Base.cross

import Base.+, Base.-, Base.*, Base./, Base.isnan, Base.inv
import Base.getindex, Base.min, Base.max
import Base.convert, Base.isapprox
import Base.size, Base.length, Base.iterate
#import Base.AbstractArray

export VectorLike, Vector3, Normal3, Point3
export Transformation, translation, scaling, rotation, look_at, swaps_handedness
export X_AXIS, Y_AXIS, Z_AXIS

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

Vector3(x::Real) = Vector3(x,x,x)
Point3(x::Real) = Point3(x,x,x)
Normal3(x::Real) = Normal3(x,x,x)

Normal3(v::Vector3) = Normal3(v.x, v.y, v.z)
Point3(v::Vector3) = Point3(v.x, v.y, v.z)

size(A::VectorLike) = 3
length(A::VectorLike) = 3

mean(A::Point3, B::Point3) = 0.5 * Point3(A.x+B.x, A.y+B.y, A.z+B.z)
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
    M::Array{Float64,2}
    MInv::Array{Float64,2}
end

const EYE = Matrix{Float64}(1.0I, 4, 4)
Transformation() = Transformation(EYE, EYE)
Transformation(M::Array{Float64,2}) = Transformation(M, inv(M))
Transformation(M::Array{T,2}) where {T<:Real} = Transformation(convert(Array{Float64,2}, M))
inv(T::Transformation) = Transformation(T.MInv, T.M)
*(T::Transformation, U::Transformation) = Transformation(T.M*U.M, U.MInv*T.MInv)
# TODO: check maths above

swaps_handedness(T::Transformation) = det(view(T.M, 1:3, 1:3)) < 0.0

translation(x::Real, y::Real, z::Real) = translation(Vector3(x,y,z))
function translation(delta::Vector3)
    M = Matrix{Float64}(1.0I, 4, 4)
    M[1,4] = delta.x
    M[2,4] = delta.y
    M[3,4] = delta.z
    return Transformation(M)
end

scaling(x::Real, y::Real, z::Real) = scaling(Vector3(x,y,z))
function scaling(scale::Vector3)
    M = Matrix{Float64}(1.0I, 4, 4)
    M[1,1] = scale.x
    M[2,2] = scale.y
    M[3,3] = scale.z
    return Transformation(M)
end

function rotation(axis::Vector3, angle::Number)
    a = axis / norm(axis)
    s = sin(angle)
    c = cos(angle)
    M = zeros(Float64, 4, 4)
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

function look_at(camera::Vector3, target::Vector3, up::Vector3)
    M = eye(4)
    dir = normalize(dir)
    left = normalize(cross(normalize(up), dir))
    newUp = normalize(cross(dir, left))
    M[1,1] = left.x 
    M[2,1] = left.y
    M[3,1] = left.z
    M[1,2] = newUp.x
    M[2,2] = newUp.y
    M[3,2] = newUp.z
    M[1,3] = dir.x
    M[2,3] = dir.y
    M[3,3] = dir.z
    M[1,4] = target.x
    M[2,4] = target.y
    M[3,4] = target.z
    return M
end

function (T::Transformation)(v::Vector3)
    return Vector3(
        T.M[1,1]*v.x + T.M[1,2]*v.y + T.M[1,3]*v.z,
        T.M[2,1]*v.x + T.M[2,2]*v.y + T.M[2,3]*v.z,
        T.M[3,1]*v.x + T.M[3,2]*v.y + T.M[3,3]*v.z
    )
end

function (T::Transformation)(v::Point3)
    invw = 1.0 / (T.M[4,1]*v.x + T.M[4,2]*v.y + T.M[4,3]*v.z + T.M[4,4])
    p = Point3((T.M[1,1]*v.x + T.M[1,2]*v.y + T.M[1,3]*v.z + T.M[1,4]) * invw,
               (T.M[2,1]*v.x + T.M[2,2]*v.y + T.M[2,3]*v.z + T.M[2,4]) * invw,
               (T.M[3,1]*v.x + T.M[3,2]*v.y + T.M[3,3]*v.z + T.M[3,4]) * invw)
    return p
end

function (T::Transformation)(v::Normal3)
    return Normal3(
        T.MInv[1,1]*v.x + T.MInv[2,1]*v.y + T.MInv[3,1]*v.z,
        T.MInv[1,2]*v.x + T.MInv[2,2]*v.y + T.MInv[3,2]*v.z,
        T.MInv[1,3]*v.x + T.MInv[2,3]*v.y + T.MInv[3,3]*v.z
    )
end

end # module
