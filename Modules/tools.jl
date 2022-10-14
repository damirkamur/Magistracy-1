module Tools

export L1_norm, L2_norm

L1_norm(vec::Vector{<:Real})::Float64 = sum((x->abs(x)).(vec))
L2_norm(vec::Vector{<:Real})::Float64 = sqrt(sum((x->x*x).(vec)))

end # module Utils