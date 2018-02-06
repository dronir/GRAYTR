
struct BoxFilter <: Filter
    xwidth::Float64
    ywidth::Float64
end

evaluate(f::BoxFilter, x, y) = 1.0


struct TriangleFilter <: Filter
    xwidth::Float64
    ywidth::Float64
end

evaluate(f::TriangleFilter, x, y) = max(0.0, (f.xwidth - abs(x)) * (f.ywidth - abs(y)))

