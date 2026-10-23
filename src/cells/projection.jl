# The projection of the Garside complex onto the cells: the splitting of the chains whose
# parent is a cell (the matching), the memoized projection p and the boundary of the cells
# (J1.d, task d3). Names provisional: the theory of the normalising homotopies of the
# article will fix the point of view.

"""
    BimoduleTerm(cell, left, right, coef)

A term `coef · left · [cell] · right` of a free bimodule over the monoid algebra of a
germ: the number of a cell (or of a chain), two simples (`0`, the identity, when there is
no multiplier) and an integer coefficient.
"""
struct BimoduleTerm
    cell::Int32
    left::Int32
    right::Int32
    coef::Int64
end

_term_key(t::BimoduleTerm) = (UInt128(UInt32(t.cell)) << 64) | (UInt128(UInt32(t.left)) << 32) | UInt128(UInt32(t.right))

# Normal form in place: sorted by (cell, left, right), equal terms summed (with a check of
# overflow), zero terms removed.
function _normalise!(buf::Vector{BimoduleTerm})
    length(buf) <= 1 && return filter!(t -> !iszero(t.coef), buf)
    sort!(buf; by = _term_key)
    k = 0
    for i in eachindex(buf)
        t = buf[i]
        if k > 0 && _term_key(buf[k]) == _term_key(t)
            s = buf[k]
            buf[k] = BimoduleTerm(s.cell, s.left, s.right, Base.checked_add(s.coef, t.coef))
        else
            k += 1
            buf[k] = t
        end
    end
    resize!(buf, k)
    return filter!(t -> !iszero(t.coef), buf)
end

# The product of two simples of G, which must be a simple (0 is the identity).
@inline function _mul(G::AbstractGerm, s::Int32, t::Int32)
    iszero(s) && return t
    iszero(t) && return s
    p = product(G, s, t)
    p === nothing && throw(ArgumentError("LCMProjectionData: the product of the simples $s and $t is not simple"))
    return Int32(p)
end

"""
    lcm_splitting(T::CellTree, h, σ) -> (y, z)

The splitting of the last letter of the chain `h|σ`, for a cell `h` of weight `ω` and a
simple `σ ≠ 1` with `ωσ` simple: `σ = y·z`, where `y` is the letter of the first child of
`h`, in the order given by the ordering, whose atom divides `ωσ` (so that `ωy = ω ∨ a`), and
`(0, σ)` (the identity for `y`) if no child of `h` has an atom dividing `ωσ`.

It classifies the chain `x = h|σ`: if `y = σ`, `x` is the corresponding child of `h`, a
cell; if `y` is the identity, `x` is collapsible and its projection is zero; otherwise `x`
is redundant, matched with the chain `h|y|z` of degree one more. The `filtfact` of the
preliminary implementation.
"""
function lcm_splitting(T::CellTree, h::Integer, σ::Integer)
    G = T.germ
    ω = cell_weight(T, h)
    ωσ = iszero(ω) ? Int32(σ) : product(G, ω, σ)
    ωσ === nothing && throw(ArgumentError("lcm_splitting: the product of the weight of the cell $h and $σ is not simple"))
    D = left_descents(G, ωσ)
    for c in cell_children(T, h)
        if T.atom[c] in D
            y = T.letter[c]
            z = ldiv(G, y, σ)
            z === nothing && throw(ArgumentError("lcm_splitting: $y does not divide $σ"))
            return (Int32(y), Int32(z))
        end
    end
    return (Int32(0), Int32(σ))
end

"""
    LCMProjectionData(T::CellTree)

The data of the projection `p` of the Garside complex onto the cells of `T`, which
[`lcm_projection!`](@ref) computes and memoizes: the tree of the cells, the tree of the
chains visited ([`ChainTree`](@ref)), and the memo, the normal form of `p(x)` for each
chain `x` visited, stored as an interval of a single arena of [`BimoduleTerm`](@ref)s.
The boundary of the cells, [`lcm_boundary`](@ref), is read from it.

The projection follows the matching of [`lcm_splitting`](@ref): `p(x) = x` for a cell;
for `x = h|σ` with `h` a cell, `p(x) = 0` if `x` is collapsible and, if `x` is redundant,
matched with `ϕx = h|y|z` of degree `n + 1`, `p(x) = (−1)ⁿ⁺¹ Σᵢ (−1)ⁱ σᵢ p(dᵢϕx) τᵢ`, the
sum over the faces of `ϕx` but the `n`-th, which is `x`; for `x = h|σ` with `h` not a
cell, `p(x) = Σ λ s p(c | t·σ)` over the terms `λ s [c] t` of `p(h)`.
"""
struct LCMProjectionData{G <: AbstractGerm, O <: LCMOrdering}
    tree::CellTree{G, O}
    chains::ChainTree
    memo_start::Vector{Int64}           # first term of p(x) in the arena
    memo_length::Vector{Int32}          # number of terms of p(x); -1 not computed, -2 in progress
    arena::Vector{BimoduleTerm}
    buffers::Vector{Vector{BimoduleTerm}}   # scratch, one for each depth of recursion
end

LCMProjectionData(T::CellTree) =
    LCMProjectionData(T, ChainTree(T), Int64[], Int32[], BimoduleTerm[], Vector{BimoduleTerm}[])

# Memo entries for the chains added to the tree since the last call.
function _sync!(P::LCMProjectionData)
    n = _nchains(P.chains)
    m = length(P.memo_length)
    if m < n
        resize!(P.memo_start, n)
        resize!(P.memo_length, n)
        for x in (m + 1):n
            P.memo_start[x] = 0
            P.memo_length[x] = -1
        end
    end
    return P
end

function _buffer!(P::LCMProjectionData, depth::Int)
    while length(P.buffers) < depth
        push!(P.buffers, BimoduleTerm[])
    end
    b = P.buffers[depth]
    empty!(b)
    return b
end

_memo_range(P::LCMProjectionData, x::Integer) =
    P.memo_start[x]:(P.memo_start[x] + P.memo_length[x] - 1)

# Appends λ · a · p(y) · b to buf, p(y) being in the memo.
function _append_scaled!(buf, P::LCMProjectionData, y::Integer, λ::Int64, a::Int32, b::Int32)
    G = P.tree.germ
    for k in _memo_range(P, y)
        t = P.arena[k]
        push!(buf, BimoduleTerm(t.cell, _mul(G, a, t.left), _mul(G, t.right, b), Base.checked_mul(λ, t.coef)))
    end
    return buf
end

# The letter k of ϕ = [h₁, …, hₙ₋₁, y, z] (hl the letters of h).
@inline _phi(hl, y::Int32, z::Int32, k::Int) = k <= length(hl) ? hl[k] : (k == length(hl) + 1 ? y : z)

"""
    lcm_projection!(P::LCMProjectionData, x) -> AbstractVector{BimoduleTerm}

The projection `p(x)` of the chain `x` of the Garside complex onto the cells, in normal
form (sorted, equal terms summed, no zero): its terms are cells of the degree of `x`, with
their left and right simples. The chain `x` is a node of `P.chains` or the list of its
letters. The result is a view of the memo, valid until the next call; the memo and the
tree of the chains grow on the way.
"""
function lcm_projection!(P::LCMProjectionData, x::Integer)
    _project!(P, Int32(x), 1)
    return view(P.arena, _memo_range(P, x))
end

function lcm_projection!(P::LCMProjectionData, letters::AbstractVector{<:Integer})
    x = _extend!(P.chains, 1, letters, eachindex(letters))
    return lcm_projection!(P, x)
end

function _project!(P::LCMProjectionData, x::Int32, depth::Int)
    _sync!(P)
    ℓ = P.memo_length[x]
    ℓ >= 0 && return nothing
    ℓ == -2 && throw(ArgumentError("LCMProjectionData: the projection of the chain $x depends on itself (the matching is not acyclic)"))
    P.memo_length[x] = -2
    C = P.chains
    T = P.tree
    G = T.germ
    buf = _buffer!(P, depth)
    if x <= C.ncells
        push!(buf, BimoduleTerm(x, 0, 0, 1))
    else
        h = C.parent[x]
        σ = C.letter[x]
        if h <= C.ncells
            y, z = lcm_splitting(T, h, σ)
            iszero(y) || _project_redundant!(buf, P, h, y, z, Int(C.degree[x]), depth)
        else
            _project!(P, h, depth + 1)
            for k in _memo_range(P, h)
                t = P.arena[k]
                c = _child!(C, t.cell, _mul(G, t.right, σ))
                _project!(P, c, depth + 1)
                _sync!(P)
                _append_scaled!(buf, P, c, t.coef, t.left, Int32(0))
            end
        end
        _normalise!(buf)
    end
    _sync!(P)
    P.memo_start[x] = length(P.arena) + 1
    P.memo_length[x] = length(buf)
    append!(P.arena, buf)
    return nothing
end

# p(x) for x = h|σ redundant, matched with ϕx = h|y|z of degree n + 1: the faces of ϕx but
# the n-th, with the sign (−1)ⁿ⁺¹.
function _project_redundant!(buf, P::LCMProjectionData, h::Int32, y::Int32, z::Int32, n::Int, depth::Int)
    C = P.chains
    G = P.tree.germ
    hl = _letters(C, h)                 # n − 1 letters
    pre = _prefixes(C, h)               # pre[k + 1]: prefix of h of length k
    sign = isodd(n + 1) ? Int64(-1) : Int64(1)
    for i in 0:(n + 1)
        i == n && continue
        λ = isodd(i) ? -sign : sign
        a = Int32(0)
        b = Int32(0)
        if i == 0                       # ϕ₁ · [ϕ₂|…|ϕₙ₊₁]
            a = _phi(hl, y, z, 1)
            node = Int32(1)
            for k in 2:(n + 1)
                node = _child!(C, node, _phi(hl, y, z, k))
            end
        elseif i == n + 1               # [ϕ₁|…|ϕₙ] · ϕₙ₊₁, that is the child h|y, times z
            node = _child!(C, h, y)
            b = z
        else                            # 1 ≤ i ≤ n − 1: [ϕ₁|…|ϕᵢϕᵢ₊₁|…|ϕₙ₊₁]
            node = _child!(C, pre[i], _mul(G, _phi(hl, y, z, i), _phi(hl, y, z, i + 1)))
            for k in (i + 2):(n + 1)
                node = _child!(C, node, _phi(hl, y, z, k))
            end
        end
        _project!(P, node, depth + 1)
        _sync!(P)
        _append_scaled!(buf, P, node, λ, a, b)
    end
    return buf
end

"""
    lcm_boundary(P::LCMProjectionData, x) -> Vector{BimoduleTerm}

The boundary of the cell `x` of degree `n ≥ 1` in the resolution: the projection of its
boundary in the Garside complex, `∂x = Σᵢ (−1)ⁱ σᵢ p(dᵢx) τᵢ`, a combination of cells of
degree `n − 1` in normal form; empty for the cell of degree 0. The column of `x` in the
matrix of the differential.
"""
lcm_boundary(P::LCMProjectionData, x::Integer) = _boundary!(BimoduleTerm[], P, x)

function _boundary!(out::Vector{BimoduleTerm}, P::LCMProjectionData, x::Integer)
    empty!(out)
    T = P.tree
    C = P.chains
    G = T.germ
    1 <= x <= ncells(T) || throw(ArgumentError("lcm_boundary: $x is not a cell"))
    n = cell_degree(T, x)
    n == 0 && return out
    xl = _letters(C, x)
    pre = _prefixes(C, x)
    for i in 0:n
        λ = isodd(i) ? Int64(-1) : Int64(1)
        a = Int32(0)
        b = Int32(0)
        if i == 0                       # x₁ · [x₂|…|xₙ]
            a = xl[1]
            node = _extend!(C, 1, xl, 2:n)
        elseif i == n                   # [x₁|…|xₙ₋₁] · xₙ, the parent of x
            node = C.parent[x]
            b = xl[n]
        else                            # [x₁|…|xᵢxᵢ₊₁|…|xₙ]
            node = _child!(C, pre[i], _mul(G, xl[i], xl[i + 1]))
            node = _extend!(C, node, xl, (i + 2):n)
        end
        _project!(P, Int32(node), 1)
        _sync!(P)
        _append_scaled!(out, P, node, λ, a, b)
    end
    return _normalise!(out)
end

function Base.show(io::IO, P::LCMProjectionData)
    print(io, "LCMProjectionData: ", ncells(P.tree), " cells, ", _nchains(P.chains),
          " chains, ", length(P.arena), " terms in the memo")
end
