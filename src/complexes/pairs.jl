# Collapsible pairs of a complex (found only; the contraction itself comes later).

const _Pair = NamedTuple{(:n, :source, :target), NTuple{3, Int}}

_is_integer_unit(c::Integer) = c == 1 || c == -1

"""
    collapsible_pairs(C::ChainComplex)

The collapsible pairs of `C`: the triples `(n, source, target)` such that the coefficient
of the `target`-th cell of `C_n` in the image of the `source`-th cell of `C_{n+1}` by
`d_n` is a unit, listed by increasing `n`, then `source`, then `target`. The number of
pairs is the length of the list, the numbers of distinct sources and targets those of the
distinct pairs `(n, source)` and `(n, target)`.

- Over a monoid algebra (coefficients in ℤ): the coefficient must be a single monomial,
  with the label of the identity and the coefficient ±1. A monomial ±γ is a unit only if
  γ is invertible, that is the identity in a Garside monoid; coefficients with several
  monomials are excluded, which loses nothing for homogeneous resolutions, where all the
  monomials of a coefficient have the same total length.
- Over ℤ: the entry must be ±1. After specialisation this finds many more pairs than in
  the resolution itself: it is a different notion, useful for the reduction of integer
  matrices rather than of resolutions.
"""
function collapsible_pairs(C::ChainComplex{A, M, LabelledSparseMatrix{L, Tv}}) where
        {A <: MonoidAlgebra, M, L, Tv <: Integer}
    base_ring(C.ring) isa AbstractAlgebra.Integers ||
        throw(ArgumentError("collapsible pairs need integer coefficients for now"))
    one_label = identity_label(L)
    pairs = _Pair[]
    for n in first(range(C)):(last(range(C)) - 1)
        d = differential(C, n)
        for j in 1:d.n
            p, b = d.colptr[j], d.colptr[j + 1] - 1
            while p <= b
                q = p                               # the monomials of row rowval[p]: p to q
                while q < b && d.rowval[q + 1] == d.rowval[p]
                    q += 1
                end
                if q == p && d.label[p] == one_label && _is_integer_unit(d.nzval[p])
                    push!(pairs, (n = n, source = j, target = Int(d.rowval[p])))
                end
                p = q + 1
            end
        end
    end
    return pairs
end

function collapsible_pairs(C::ChainComplex{R, M, D}) where
        {R <: AbstractAlgebra.Integers, M, D <: SparseMatrixCSC{<:Integer}}
    pairs = _Pair[]
    for n in first(range(C)):(last(range(C)) - 1)
        d = differential(C, n)
        rows, vals = rowvals(d), nonzeros(d)
        for j in 1:size(d, 2), k in nzrange(d, j)
            _is_integer_unit(vals[k]) && push!(pairs, (n = n, source = j, target = Int(rows[k])))
        end
    end
    return pairs
end
