# Checked arithmetic on coefficients: machine integers raise an OverflowError instead of
# wrapping around silently; other types (BigInt, ...) use their own operations.

_checked_add(a::T, b::T) where {T <: Base.BitInteger} = Base.checked_add(a, b)
_checked_add(a, b) = a + b

_checked_mul(a::T, b::T) where {T <: Base.BitInteger} = Base.checked_mul(a, b)
_checked_mul(a, b) = a * b
