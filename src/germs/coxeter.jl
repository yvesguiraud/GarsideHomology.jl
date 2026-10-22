# Lazy germs of the classical braid monoids B⁺(W) of the finite crystallographic Coxeter
# groups, computed on the root system.

"""
    CoxeterGerm{n} <: AbstractGerm
    CoxeterGerm(C)

The germ of the classical braid monoid B⁺(W) of the finite Coxeter group W of rank `n`
given by its integral Cartan matrix `C` (crystallographic types), the simples being the
elements of W. The germ is lazy: its simples are numbered as they are met (identity
`0`, atoms `1:n` in the order of the rows of `C`, the others in the order of the queries,
which is not canonical), for the large groups that cannot be tabulated, such as E₈
(696 729 600 simples, of which the construction of the resolution meets a few percent).

An element w is identified by its key, the indices of the roots w(α₁), …, w(αₙ), with
the convention of Chevie for the product (`σ·τ` acts as σ then τ on the indices of the
roots). For each simple met, only its key and its length are kept; left descents are
read on the key, and the other operations computed with the coordinates of the roots
in the basis of the simple roots. The lcms with the atoms are memorised. A query numbers
only the simple it returns: divisions, words and lexicographic comparisons follow the
smallest words on the keys, without numbering the intermediate quotients.

[`nsimples`](@ref) is |W|, from the exponents (Kostant: heights of the positive roots);
[`nnumbered`](@ref) is the number of simples numbered so far; [`simples`](@ref) is not
defined. A query may number new simples: a `CoxeterGerm` must not be queried from
several threads at once. Two `CoxeterGerm` are equal only if identical, their
numberings depending on the history of their queries. `C` and its transpose give the
same germ (the Weyl groups of dual root systems are the same Coxeter group).
"""
struct CoxeterGerm{n} <: AbstractGerm
    cartan::Matrix{Int}
    npositive::Int                           # N; roots 1:N positive, N + k = −(root k)
    coords::Vector{NTuple{n, Int8}}          # coordinates of the roots, simple roots 1:n
    heights::Vector{Int16}                   # heights of the roots (negative for negative roots)
    root_index::Dict{NTuple{n, Int8}, Int16}
    reflect::Vector{Vector{Int16}}           # reflect[i][r]: index of sᵢ(root r)
    coxeter::Matrix{Int}                     # Coxeter matrix m_ab
    order::Int                               # |W|
    keys::Vector{NTuple{n, Int16}}           # keys[s + 1]
    lengths::Vector{UInt8}                   # lengths[s + 1]
    index::Dict{NTuple{n, Int16}, Int32}     # key → simple
    lcm_memo::Dict{Tuple{Int32, Int32}, NTuple{2, Int32}}
end

"""
    nnumbered(G)

The number of simples numbered so far by the lazy germ `G`; it only grows.
"""
function nnumbered end

function CoxeterGerm(C::AbstractMatrix{<:Integer})
    n = size(C, 1)
    size(C) == (n, n) && n >= 1 || throw(ArgumentError("CoxeterGerm: the Cartan matrix must be square"))
    n <= 64 || throw(ArgumentError("CoxeterGerm: rank $n larger than 64"))
    return _coxeter_germ(Val(n), Matrix{Int}(C))
end

function _coxeter_germ(::Val{n}, C::Matrix{Int}) where {n}
    for i in 1:n, j in 1:n
        if i == j
            C[i, i] == 2 || throw(ArgumentError("CoxeterGerm: diagonal entry C[$i, $i] = $(C[i, i]) ≠ 2"))
        else
            C[i, j] <= 0 && (C[i, j] == 0) == (C[j, i] == 0) ||
                throw(ArgumentError("CoxeterGerm: not a Cartan matrix at ($i, $j)"))
        end
    end
    m = [i == j ? 1 : _coxeter_order(C[i, j] * C[j, i]) for i in 1:n, j in 1:n]
    # Positive roots by closure of the simple roots under the simple reflections:
    # sᵢ(β) = β − ⟨β, αᵢ^∨⟩ αᵢ with ⟨β, αᵢ^∨⟩ = Σⱼ C[i, j] βⱼ.
    simple = [ntuple(k -> k == i ? 1 : 0, n) for i in 1:n]
    seen = Set(simple)
    queue = copy(simple)
    head = 1
    while head <= length(queue)
        β = queue[head]
        head += 1
        for i in 1:n
            β == simple[i] && continue
            c = sum(C[i, j] * β[j] for j in 1:n)
            γ = ntuple(k -> k == i ? β[k] - c : β[k], n)
            all(>=(0), γ) && !(γ in seen) || continue
            push!(seen, γ)
            push!(queue, γ)
            length(queue) <= 10_000 ||
                throw(ArgumentError("CoxeterGerm: more than 10 000 positive roots: not of finite type"))
        end
    end
    others = sort!(setdiff(queue, simple); by = β -> (sum(β), β))
    positive = vcat(simple, others)
    N = length(positive)
    maximum(maximum(abs, β) for β in positive) <= typemax(Int8) ||
        throw(ArgumentError("CoxeterGerm: root coordinates too large"))
    coords = [ntuple(k -> Int8(β[k]), n) for β in vcat(positive, [map(-, β) for β in positive])]
    heights = Int16[sum(Int, c) for c in coords]
    root_index = Dict{NTuple{n, Int8}, Int16}(c => Int16(r) for (r, c) in enumerate(coords))
    reflect = [Int16[_root_of(root_index, _reflection(C, i, coords[r], Val(n))) for r in 1:(2N)] for i in 1:n]
    order = _weyl_order(C, positive)
    identity_key = ntuple(j -> Int16(j), Val(n))
    G = CoxeterGerm{n}(C, N, coords, heights, root_index, reflect, m, order,
                       [identity_key], UInt8[0], Dict(identity_key => Int32(0)),
                       Dict{Tuple{Int32, Int32}, NTuple{2, Int32}}())
    for i in 1:n
        _register!(G, ntuple(j -> reflect[i][j], Val(n)), 1) == i ||
            throw(ArgumentError("CoxeterGerm: two simple reflections coincide"))
    end
    return G
end

function _coxeter_order(p)
    p == 0 && return 2
    p == 1 && return 3
    p == 2 && return 4
    p == 3 && return 6
    throw(ArgumentError("CoxeterGerm: Cᵢⱼ·Cⱼᵢ = $p is not 0, 1, 2 or 3 (not of finite crystallographic type)"))
end

_reflection(C, i, β::NTuple{n}, ::Val{n}) where {n} =
    (c = sum(C[i, j] * β[j] for j in 1:n); ntuple(k -> k == i ? β[k] - c : β[k], Val(n)))

function _root_of(root_index, v)
    r = get(root_index, map(Int8, v), Int16(0))
    r == 0 && throw(ArgumentError("CoxeterGerm: $v is not a root"))
    return r
end

# |W|, from the exponents of each irreducible component: the dual partition of the
# numbers of positive roots of each height (Kostant), degrees = exponents + 1.
function _weyl_order(C, positive)
    n = size(C, 1)
    component = collect(1:n)
    root(i) = component[i] == i ? i : (component[i] = root(component[i]))
    for i in 1:n, j in 1:n
        i != j && C[i, j] != 0 && (component[root(i)] = root(j))
    end
    order = big(1)
    for r in unique(root(i) for i in 1:n)
        heights = [sum(β) for β in positive if root(findfirst(!iszero, β)) == r]
        counts = [count(==(k), heights) for k in 1:maximum(heights)]
        rank = counts[1]
        exponents = [count(>=(j), counts) for j in 1:rank]
        order *= prod(big.(exponents) .+ 1)
    end
    order <= typemax(Int) || throw(ArgumentError("CoxeterGerm: |W| = $order too large"))
    return Int(order)
end

# ------------------------------------------------------------------------------------
# Numbering

function _register!(G::CoxeterGerm{n}, k::NTuple{n, Int16}, ℓ::Integer) where {n}
    id = get(G.index, k, Int32(-1))
    id >= 0 && return id
    m = length(G.keys)
    m < typemax(Int32) || throw(ArgumentError("CoxeterGerm: more than typemax(Int32) simples"))
    ℓ <= 255 || throw(ArgumentError("CoxeterGerm: length $ℓ larger than 255"))
    push!(G.keys, k)
    push!(G.lengths, UInt8(ℓ))
    G.index[k] = Int32(m)
    return Int32(m)
end

_key(G::CoxeterGerm, s::Integer) = G.keys[Int(s) + 1]

# The key of σ·a, a an atom: (j^σ)^a.
_right_atom(G::CoxeterGerm, k::NTuple, a::Integer) = map(r -> G.reflect[a][r], k)

# The key of a·σ, a an atom: the root σ(αⱼ − ⟨αⱼ, αₐ^∨⟩ αₐ).
function _left_atom(G::CoxeterGerm{n}, k::NTuple{n, Int16}, a::Integer) where {n}
    ca = G.coords[k[a]]
    return ntuple(Val(n)) do j
        c = G.cartan[a, j]
        cj = G.coords[k[j]]
        G.root_index[ntuple(i -> Int8(cj[i] - c * ca[i]), Val(n))]
    end
end

# The key of σ·τ: (j^σ)^τ = τ applied to the root j^σ, by linearity.
function _product_key(G::CoxeterGerm{n}, kσ::NTuple{n, Int16}, kτ::NTuple{n, Int16}) where {n}
    return ntuple(Val(n)) do j
        c = G.coords[kσ[j]]
        v = ntuple(i -> sum(Int(c[l]) * G.coords[kτ[l]][i] for l in 1:n), Val(n))
        G.root_index[map(Int8, v)]
    end
end

# The length of the element of key k: the positive roots sent to negative roots.
function _length_of(G::CoxeterGerm{n}, k::NTuple{n, Int16}) where {n}
    h = map(r -> Int(G.heights[r]), k)
    ℓ = 0
    for r in 1:G.npositive
        β = G.coords[r]
        sum(Int(β[l]) * h[l] for l in 1:n) < 0 && (ℓ += 1)
    end
    return ℓ
end

# ------------------------------------------------------------------------------------
# Interface

natoms(G::CoxeterGerm{n}) where {n} = n
nsimples(G::CoxeterGerm) = G.order
nnumbered(G::CoxeterGerm) = length(G.keys)
Base.length(G::CoxeterGerm, s::Integer) = Int(G.lengths[Int(s) + 1])

function left_descents(G::CoxeterGerm{n}, s::Integer) where {n}
    k = _key(G, s)
    N = G.npositive
    return StaticBitSet{1}(i for i in 1:n if k[i] > N)
end

function product(G::CoxeterGerm{n}, s::Integer, t::Integer) where {n}
    x, y = Int32(s), Int32(t)
    iszero(y) && (_key(G, x); return x)
    iszero(x) && (_key(G, y); return y)
    ℓ = length(G, x) + length(G, y)
    k = y <= n ? _right_atom(G, _key(G, x), y) : _product_key(G, _key(G, x), _key(G, y))
    id = get(G.index, k, Int32(-1))
    if id >= 0
        return length(G, id) == ℓ ? id : nothing
    end
    return _length_of(G, k) == ℓ ? _register!(G, k, ℓ) : nothing
end

# The smallest left descent of the element of key k, 0 for the identity (a loop rather than
# findfirst with a closure, which would box a key reassigned in the callers' loops).
@inline function _first_descent(k::NTuple{n, Int16}, N) where {n}
    for i in 1:n
        k[i] > N && return i
    end
    return 0
end

# s\t along the smallest word of s, on the keys: only the result is numbered, not the
# intermediate quotients a\s and a\t.
function ldiv(G::CoxeterGerm{n}, s::Integer, t::Integer) where {n}
    x, y = Int32(s), Int32(t)
    kx = _key(G, x)
    iszero(x) && (_key(G, y); return y)
    x == y && return Int32(0)
    ℓx, ℓy = length(G, x), length(G, y)
    ℓx < ℓy || return nothing
    ky = _key(G, y)
    N = G.npositive
    for _ in 1:ℓx
        a = _first_descent(kx, N)                          # smallest left descent of x
        ky[a] > N || return nothing                        # a ⋠ y
        kx = _left_atom(G, kx, a)
        ky = _left_atom(G, ky, a)
    end
    return _register!(G, ky, ℓy - ℓx)
end

# The smallest word and the lexicographic comparison, on the keys, numbering nothing.
function word(G::CoxeterGerm{n}, s::Integer) where {n}
    k = _key(G, s)
    N = G.npositive
    w = Vector{Int32}(undef, length(G, s))
    for i in eachindex(w)
        a = _first_descent(k, N)
        w[i] = a
        k = _left_atom(G, k, a)
    end
    return w
end

function isless_lex(G::CoxeterGerm{n}, s::Integer, t::Integer) where {n}
    kx, ky = _key(G, s), _key(G, t)
    e = G.keys[1]
    N = G.npositive
    while true
        kx == ky && return false
        ky == e && return false
        kx == e && return true
        a = _first_descent(kx, N)
        b = _first_descent(ky, N)
        a == b || return a < b
        kx = _left_atom(G, kx, a)
        ky = _left_atom(G, ky, b)
    end
end

function lcm_atom(G::CoxeterGerm{n}, ω::Integer, a::Integer) where {n}
    1 <= a <= n || throw(BoundsError(1:n, a))
    x = Int32(ω)
    _key(G, x)
    m, τ = get!(G.lcm_memo, (x, Int32(a))) do
        m = _lcm(G, x, Int32(a))
        (m, ldiv(G, x, m)::Int32)
    end
    return (lcm = m, complement = τ)
end

# Right lcm of the simples x and y, by the recursion of Garside theory, from the lcms of
# two atoms: a ∨ b = a·(a b a … of length m_ab − 1 starting with b).
function _lcm(G::CoxeterGerm{n}, x::Int32, y::Int32) where {n}
    is_left_divisor(G, y, x) && return x
    is_left_divisor(G, x, y) && return y
    if y <= n                                           # an atom a = y, a ⋠ x
        a = y
        b = Int32(first(left_descents(G, x)))           # x = b·x′, b ≠ a
        x′ = ldiv(G, b, x)::Int32
        c = _alternating(G, a, b, G.coxeter[a, b] - 1)  # b ∨ a = b·c
        return product(G, b, _lcm(G, x′, c))::Int32
    end
    y₁ = Int32(first(left_descents(G, y)))             # y = y₁·y′
    y′ = ldiv(G, y₁, y)::Int32
    z = _lcm(G, x, y₁)                                  # x ∨ y = (x ∨ y₁) ∨ y
    z′ = ldiv(G, y₁, z)::Int32
    return product(G, y₁, _lcm(G, z′, y′))::Int32
end

# The simple a b a … with `len` letters.
function _alternating(G::CoxeterGerm, a::Int32, b::Int32, len::Integer)
    s = Int32(0)
    for i in 1:len
        s = product(G, s, isodd(i) ? a : b)::Int32
    end
    return s
end

opposite(G::CoxeterGerm) = CoxeterGerm(G.cartan)

Base.show(io::IO, G::CoxeterGerm{n}) where {n} =
    print(io, "CoxeterGerm of rank ", n, ", ", nnumbered(G), " of ", nsimples(G), " simples numbered")
