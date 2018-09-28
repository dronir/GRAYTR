
"""
    lerp(t, a, b)

Linear interpolation on the interval `[a, b]` as `t ∈ [0, 1].
"""
lerp(t, a, b) = (one(t)-t)*a + t*b


"""
    quadratic(A, B, C)

Solve quadratic equation `A x^2 + B x + C = 0` for real coefficients.

If a real solution exists, returns `(true, t0, t1)` where `t0` is the smaller
and `t1` the larger solution. If no real solution exists, returns `(false, Inf, Inf)`.
"""
function quadratic(A::Real, B::Real, C::Real)
    dd = B^2 - 4*A*C
    if dd < 0.0
        return false, Inf, Inf
    end
    q = B<0 ? -0.5 * (B - sqrt(dd)) : -0.5 * (B + sqrt(dd))
    t0 = q/A
    t1 = C/q
    return true, min(t0,t1), max(t0,t1)
end


"""
    normal_derivatives(dpdu, dpdv, dpduu, dpduv, dpdvv)

Compute partial derivatives of normal vectors from Weingarten equations,
given first and second derivatives of position vector.
Used by several shape models.
"""
function normal_derivatives(dpdu, dpdv, dpduu, dpduv, dpdvv)
    E = dot(dpdu, dpdu)
    F = dot(dpdu, dpdv)
    G = dot(dpdv, dpdv)
    N = normalize(cross(dpdu, dpdv))
    e = dot(N, dpduu)
    f = dot(N, dpduv)
    g = dot(N, dpdvv)
    
    invEG = 1.0 / (E*G - F^2)
    dndu = Normal3((f*F - e*G) * invEG * dpdu + (e*F - f*E) * invEG * dpdv)
    dndv = Normal3((g*F - f*G) * invEG * dpdu + (f*F - g*E) * invEG * dpdv)
    return dndu, dndv
end

