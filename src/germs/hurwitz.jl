# The Hurwitz action of the braid group on tuples, independent of any representation of
# the group elements.

"""
    hurwitz_orbit(start::NTuple{n, <:Integer}, conjugate; limit = 10^8) -> Vector

The orbit of the tuple `start` under the Hurwitz action of the braid group on `n`
strands: the generator σₖ replaces the entries `(x, y)` at the positions `k, k + 1` by
`(x y x⁻¹, x)`, which preserves the product of the tuple. The entries are indices of
group elements, and `conjugate(i, j)` must return the index of xᵢ xⱼ xᵢ⁻¹.

The orbit is computed by a breadth-first search with the moves σₖ only: this suffices
for a finite orbit, on which each σₖ acts as a permutation. It is returned in the order
of discovery, `start` first; an `ArgumentError` is thrown if it exceeds `limit` tuples.
For a Coxeter element c of a well-generated complex reflection group and a reduced
decomposition of c into reflections, the orbit is the set of all reduced decompositions
of c (Bessis), of size n!·hⁿ/|W|.
"""
function hurwitz_orbit(start::NTuple{n, T}, conjugate; limit::Integer = 10^8) where {n, T <: Integer}
    seen = Set{NTuple{n, T}}((start,))
    orbit = NTuple{n, T}[start]
    head = 1
    while head <= length(orbit)
        t = orbit[head]
        head += 1
        for k in 1:(n - 1)
            u = Base.setindex(Base.setindex(t, T(conjugate(t[k], t[k + 1])), k), t[k], k + 1)
            u in seen && continue
            push!(seen, u)
            push!(orbit, u)
            length(orbit) <= limit ||
                throw(ArgumentError("hurwitz_orbit: more than $limit tuples"))
        end
    end
    return orbit
end
