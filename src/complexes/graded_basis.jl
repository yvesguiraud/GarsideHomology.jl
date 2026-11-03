# The graded basis of a complex, and the printing of long lists within a line (J1.h). The
# generators are read from their written form with the elements (elements.jl).

"""
    GradedBasis

The graded basis of a chain complex `C`, as returned by `basis(C)`: `X[n]` is the basis of
`C_n`, `basis(C[n])`, for `n` in `range(X)`, the degrees of `C` (a `BoundsError` outside
them); `X[end]`, `length(X)` (the number of degrees), `pairs(X)` (degree and basis) and the
iteration over the bases, by increasing degree. At the REPL, one line per degree, cut in the
middle when it is too long; `show(io, X)` gives the summary line.
"""
struct GradedBasis{C <: AbstractChainComplex}
    complex::C
end

"""
    basis(C::AbstractChainComplex) -> GradedBasis

The graded basis of the complex `C`, a [`GradedBasis`](@ref): `basis(C)[n] == basis(C[n])`.
A method of the function `basis` of AbstractAlgebra, as `basis(M)` for a module.
"""
AbstractAlgebra.basis(C::AbstractChainComplex) = GradedBasis(C)

Base.range(X::GradedBasis) = range(X.complex)
Base.getindex(X::GradedBasis, n::Integer) = basis(X.complex[n])
Base.firstindex(X::GradedBasis) = first(range(X))
Base.lastindex(X::GradedBasis) = last(range(X))
Base.length(X::GradedBasis) = length(range(X))
Base.keys(X::GradedBasis) = range(X)

function Base.iterate(X::GradedBasis, n::Int = first(range(X)))
    n > last(range(X)) && return nothing
    return X[n], n + 1
end

Base.:(==)(X::GradedBasis, Y::GradedBasis) = range(X) == range(Y) && all(X[n] == Y[n] for n in range(X))
Base.hash(X::GradedBasis, h::UInt) = hash(range(X), hash(:GradedBasis, h))

Base.show(io::IO, X::GradedBasis) =
    print(io, "GradedBasis in degrees ", range(X), ", ranks ", [length(X[n]) for n in range(X)])

function Base.show(io::IO, ::MIME"text/plain", X::GradedBasis)
    print(io, "GradedBasis in degrees ", range(X))
    width = displaysize(io)[2]
    for n in range(X)
        prefix = "  $n: "
        print(io, "\n", prefix)
        _print_limited(io, X[n], width - textwidth(prefix))
    end
end

# Prints the elements of v separated by ", ". When io has :limit, the line keeps within
# `width` characters: as many elements as fit, taken alternately at the head and at the
# tail, with "…" between them (the first element is always written).
function _print_limited(io::IO, v::AbstractVector, width::Integer)
    n = length(v)
    text(k) = sprint(print, v[k]; context = io)
    if !(get(io, :limit, false)::Bool) || _fits(text, n, width)
        for k in 1:n
            k > 1 && print(io, ", ")
            print(io, text(k))
        end
        return nothing
    end
    head, tail = String[], String[]
    used = 0                                 # widths of the elements kept and their separators
    i, j = 1, n
    from_head = true
    while i <= j
        s = text(from_head ? i : j)
        separator = isempty(head) && isempty(tail) ? 0 : 2
        reserve = j > i ? 5 : 0              # ", …, " if elements remain after this one
        isempty(head) || used + separator + textwidth(s) + reserve <= width || break
        used += separator + textwidth(s)
        if from_head
            push!(head, s)
            i += 1
        else
            pushfirst!(tail, s)
            j -= 1
        end
        from_head = !from_head
    end
    # Not all kept, since they do not fit (only the first one when it is too wide alone).
    print(io, join(head, ", "))
    i <= j && print(io, isempty(tail) ? ", …" : ", …, ")
    print(io, join(tail, ", "))
    return nothing
end

# Whether the n elements, separated by ", ", hold in `width` characters; only the elements
# up to the first overflow are written.
function _fits(text, n::Integer, width::Integer)
    total = 0
    for k in 1:n
        total += textwidth(text(k)) + (k > 1 ? 2 : 0)
        total > width && return false
    end
    return true
end
