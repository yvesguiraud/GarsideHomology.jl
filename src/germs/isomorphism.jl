# Isomorphisms of germs: a bijection of the atoms extending to the simples.

"""
    germ_isomorphism(G::AbstractGerm, H::AbstractGerm; limit = 10^6)
        -> Union{Vector{Int32}, Nothing}

An isomorphism of germs from `G` to `H`, given by the images of the atoms: `φ[a]` is the
atom of `H` image of the atom `a` of `G`; `nothing` when the germs are not isomorphic. An
isomorphism is a bijection of the atoms that extends to a bijection of the simples with
φ(σ·a) = φ(σ)·φ(a), the product being defined in `G` exactly when it is defined in `H`;
it then preserves the lengths, the divisibility, the lcms and Δ.

The search is a backtracking on the atoms of `G`, taken by decreasing number of defined
products a·b of two atoms (the atom of smaller number first on a tie), an atom being sent only to
an atom of `H` with the same numbers of defined products on the right and on the left.
The images of the products of two atoms already placed are kept consistent along the way,
which imposes the relations ab = cd. Once all the atoms are placed, the map is extended
to the simples along the right covers, from the identity, and checked on all of them.
`limit` bounds the number of trials of an image of an atom: an `ArgumentError` is thrown
when it is reached, so that `nothing` always means that the germs are not isomorphic.

Germs with different numbers of atoms or of simples are told apart without search. Every
simple of `G` and `H` is visited by products with the atoms: for a lazy germ such as
[`CoxeterGerm`](@ref), all its simples get numbered.
"""
function germ_isomorphism(G::AbstractGerm, H::AbstractGerm; limit::Integer = 10^6)
    limit >= 1 || throw(ArgumentError("germ_isomorphism: limit = $limit < 1"))
    N = natoms(G)
    (natoms(H) == N && nsimples(G) == nsimples(H)) || return nothing
    PG, PH = _atom_products(G), _atom_products(H)
    rightG, rightH = _defined_counts(PG, 2), _defined_counts(PH, 2)
    leftG, leftH = _defined_counts(PG, 1), _defined_counts(PH, 1)
    order = sortperm(rightG; rev = true)  # stable: the smaller number first on a tie
    search = _IsomorphismSearch(PG, PH, rightG, rightH, leftG, leftH, order,
                                zeros(Int32, N), falses(N), zeros(Int32, nsimples(G)),
                                zeros(Int32, nsimples(H)), Int32[], 0, Int(limit))
    _place_atom!(search, G, H, 1) || return nothing
    return search.φ
end

# The products of two atoms: P[a, b] = a·b, or -1 when it is not simple.
function _atom_products(G::AbstractGerm)
    N = natoms(G)
    P = fill(Int32(-1), N, N)
    for a in 1:N, b in 1:N
        p = product(G, a, b)
        p === nothing || (P[a, b] = p)
    end
    return P
end

# The number of defined products of each atom: on the right (a·b, dims = 2) or on the
# left (b·a, dims = 1).
_defined_counts(P::Matrix{Int32}, dims::Integer) = vec(count(>=(0), P; dims))

# State of the search. `φ[a]` is the image of the atom a (0 if not placed), `used[b]`
# whether the atom b of H is an image; `image[σ + 1]` is the image of the simple σ of
# length 2 of G (0 if none, the identity being never one), `preimage[τ + 1]` the inverse,
# and `added` the simples of length 2 whose images were set, in order, for undoing.
mutable struct _IsomorphismSearch
    PG::Matrix{Int32}
    PH::Matrix{Int32}
    rightG::Vector{Int}
    rightH::Vector{Int}
    leftG::Vector{Int}
    leftH::Vector{Int}
    order::Vector{Int}
    φ::Vector{Int32}
    used::BitVector
    image::Vector{Int32}
    preimage::Vector{Int32}
    added::Vector{Int32}
    trials::Int
    limit::Int
end

# Places the atoms order[t:end], the previous ones being placed; true once an isomorphism
# is found (left in search.φ), false when there is none extending the current placement.
function _place_atom!(search::_IsomorphismSearch, G::AbstractGerm, H::AbstractGerm, t::Int)
    N = length(search.φ)
    t > N && return _extends_to_simples(search.φ, G, H)
    a = search.order[t]
    for b in 1:N
        (search.used[b] || search.rightH[b] != search.rightG[a] ||
         search.leftH[b] != search.leftG[a]) && continue
        search.trials += 1
        search.trials > search.limit &&
            throw(ArgumentError("germ_isomorphism: search limit $(search.limit) reached"))
        mark = length(search.added)
        if _consistent!(search, a, b)
            search.φ[a] = b
            search.used[b] = true
            _place_atom!(search, G, H, t + 1) && return true
            search.φ[a] = 0
            search.used[b] = false
        end
        _undo!(search, mark)
    end
    return false
end

# Whether sending a to b keeps the products of two placed atoms consistent (a·c and c·a
# for every placed c, and a·a): defined on both sides together, with images forming a
# bijection; the new images are recorded in search.added.
function _consistent!(search::_IsomorphismSearch, a::Int, b::Int)
    for c in 1:length(search.φ)
        d = c == a ? b : Int(search.φ[c])
        d == 0 && continue
        _match!(search, search.PG[a, c], search.PH[b, d]) || return false
        _match!(search, search.PG[c, a], search.PH[d, b]) || return false
    end
    return true
end

function _match!(search::_IsomorphismSearch, x::Int32, y::Int32)
    (x < 0) == (y < 0) || return false
    x < 0 && return true
    current = search.image[x + 1]
    current != 0 && return current == y
    search.preimage[y + 1] != 0 && return false
    search.image[x + 1] = y
    search.preimage[y + 1] = x
    push!(search.added, x)
    return true
end

function _undo!(search::_IsomorphismSearch, mark::Int)
    while length(search.added) > mark
        x = pop!(search.added)
        search.preimage[search.image[x + 1] + 1] = 0
        search.image[x + 1] = 0
    end
    return nothing
end

# Whether the bijection φ of the atoms extends to the simples: from the identity, along
# the right covers σ → σ·a of G, the image of σ·a being φ(σ)·φ(a), defined exactly when
# σ·a is, the same along every path, and never the image of another simple.
function _extends_to_simples(φ::Vector{Int32}, G::AbstractGerm, H::AbstractGerm)
    N = length(φ)
    f = fill(Int32(-1), nsimples(G))   # f[σ + 1] = φ(σ), or -1
    taken = falses(nsimples(H))
    f[one(G) + 1] = one(H)
    taken[one(H) + 1] = true
    layer, next = Int32[one(G)], Int32[]
    reached = 1
    while !isempty(layer)
        empty!(next)
        for σ in layer, a in 1:N
            p = product(G, σ, a)
            q = product(H, f[σ + 1], φ[a])
            (p === nothing) == (q === nothing) || return false
            p === nothing && continue
            if f[p + 1] >= 0
                f[p + 1] == q || return false
            else
                taken[q + 1] && return false
                f[p + 1] = q
                taken[q + 1] = true
                reached += 1
                push!(next, p)
            end
        end
        layer, next = next, layer
    end
    return reached == nsimples(G)
end
