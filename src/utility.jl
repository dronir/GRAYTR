
const LIGHT_SPEED = 299792458 # m/s^2

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
    round_pow2(n::Integer)

Round an integer up to nearest power of two.
"""
round_pow2(n::Integer) = 2^convert(typeof(n), ceil(log(2, n)))

