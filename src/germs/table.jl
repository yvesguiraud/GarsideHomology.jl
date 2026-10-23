# Tabulated germs: every simple numbered, dense tables of size N × n.

"""
    TableGerm{K} <: AbstractGerm

A germ all of whose simples are numbered in the canonical order (identity `0`, atoms
`1:N`, then by length and smallest word), with dense tables of size `N × n` (`N` atoms,
`n` simples): the right covers `sa`, the left quotients `a \\ s` by the atoms and the
lcms `s ∨ a` with their complements; the left descents (a [`StaticBitSet{K}`](@ref)),
the lengths and the lexicographic ranks of the smallest words. Products and left
divisions of any two simples are computed along the smallest word of one of them, in
`ℓ` steps, with direct answers for the identity and the atoms. The memory is linear in
`n·N`: about 16 bytes per simple and atom, 49 MB for B*(E₈).

Built by [`tabulate`](@ref) or `TableGerm(R, identity)` from the right covers, or by
`TableGerm(G)` from any germ whose simples are all numbered. The tables are written only
by the constructor, so that a `TableGerm` can be read by several threads at once.
"""
struct TableGerm{K} <: AbstractGerm
    natoms::Int
    lengths::Vector{UInt8}               # lengths[s + 1] = ℓ(s)
    descents::Vector{StaticBitSet{K}}    # descents[s + 1]
    right_covers::Matrix{Int32}          # right_covers[a, s + 1] = sa, or -1
    left_quotients::Matrix{Int32}        # left_quotients[a, ν + 1] = a \ ν, or -1
    atom_lcms::Matrix{NTuple{2, Int32}}  # atom_lcms[a, ω + 1] = (ω ∨ a, ω \ (ω ∨ a))
    lex_ranks::Vector{Int32}             # lex_ranks[s + 1]
end

@inline _ix(s::Integer) = Int(s) + 1

# ------------------------------------------------------------------------------------
# Interface

natoms(G::TableGerm) = G.natoms
nsimples(G::TableGerm) = length(G.lengths)
simples(G::TableGerm) = Int32(0):Int32(length(G.lengths) - 1)
Base.length(G::TableGerm, s::Integer) = Int(G.lengths[_ix(s)])
left_descents(G::TableGerm, s::Integer) = G.descents[_ix(s)]
isless_lex(G::TableGerm, s::Integer, t::Integer) = G.lex_ranks[_ix(s)] < G.lex_ranks[_ix(t)]
# The canonical numbering is the degree-lexicographic order (bounds checked as elsewhere).
isless_deglex(G::TableGerm, s::Integer, t::Integer) =
    (checkbounds(G.lengths, _ix(s)); checkbounds(G.lengths, _ix(t)); Int32(s) < Int32(t))

function lcm_atom(G::TableGerm, ω::Integer, a::Integer)
    m, τ = G.atom_lcms[a, _ix(ω)]
    return (lcm = m, complement = τ)
end

function product(G::TableGerm, s::Integer, t::Integer)
    checkbounds(G.lengths, _ix(s))
    checkbounds(G.lengths, _ix(t))
    x, y = Int32(s), Int32(t)
    iszero(y) && return x
    iszero(x) && return y
    R = G.right_covers
    if y <= G.natoms                              # an atom: one lookup
        v = R[y, _ix(x)]
        return v < 0 ? nothing : v
    end
    Q = G.left_quotients
    while !iszero(y)                              # along the smallest word of t
        a = first(G.descents[_ix(y)])
        x = R[a, _ix(x)]
        x < 0 && return nothing
        y = Q[a, _ix(y)]
    end
    return x
end

function ldiv(G::TableGerm, s::Integer, ν::Integer)
    checkbounds(G.lengths, _ix(s))
    checkbounds(G.lengths, _ix(ν))
    x, y = Int32(s), Int32(ν)
    iszero(x) && return y
    x == y && return Int32(0)
    Q = G.left_quotients
    if x <= G.natoms                              # an atom: one lookup
        v = Q[x, _ix(y)]
        return v < 0 ? nothing : v
    end
    G.lengths[_ix(x)] < G.lengths[_ix(y)] || return nothing
    while !iszero(x)                              # along the smallest word of s
        a = first(G.descents[_ix(x)])
        y = Q[a, _ix(y)]
        y < 0 && return nothing
        x = Q[a, _ix(x)]
    end
    return y
end

# ------------------------------------------------------------------------------------
# Equality, display

Base.:(==)(G::TableGerm, H::TableGerm) =
    G.natoms == H.natoms && G.lengths == H.lengths && G.descents == H.descents &&
    G.right_covers == H.right_covers && G.left_quotients == H.left_quotients &&
    G.atom_lcms == H.atom_lcms && G.lex_ranks == H.lex_ranks
Base.hash(G::TableGerm, h::UInt) =
    hash(G.right_covers, hash(G.lengths, hash(G.natoms, hash(:TableGerm, h))))

Base.show(io::IO, G::TableGerm) =
    print(io, "TableGerm with ", G.natoms, " atoms and ", length(G.lengths), " simples")

function Base.show(io::IO, ::MIME"text/plain", G::TableGerm{K}) where {K}
    println(io, "TableGerm{", K, "} with ", G.natoms, " atoms and ", length(G.lengths), " simples")
    print(io, "  maximal length ", Int(maximum(G.lengths)), ", ",
          round(Base.summarysize(G) / 2^20; digits = 1), " MiB")
end

# ------------------------------------------------------------------------------------
# Construction from the right covers

_tab_error(msg) = throw(ArgumentError("tabulate: " * msg))

"""
    tabulate(R, identity) -> (G, renumber)

The [`TableGerm`](@ref) `G` of the germ given by its right covers, and the numbering of
its simples. The simples are numbered `1:n` in any order (a provisional numbering),
`identity` being the identity; `R` is an `N × n` matrix of integers whose entry
`R[a, σ]` is the simple `σa` if it is simple, `0` otherwise, so that the atoms are the
`R[a, identity]`, in the order of the rows. `renumber[p]` is the number in `G` of the
simple numbered `p` in `R`; the inverse numbering is `sortperm(renumber)`.

The construction checks (and throws an `ArgumentError` naming the provisional simples
otherwise): the entries and the atoms; that every simple is reached from the identity;
homogeneity (each cover adds 1 to the length, lengths at most 255); right and left
cancellation; associativity; that the germ is bounded (a single simple without right
cover); that each simple has a common multiple with each atom, the smallest of which is
unique among the common multiples of its length. The full lattice property and the
numbers of simples and of intervals are checked separately.
"""
tabulate(R0::AbstractMatrix{<:Integer}, identity::Integer) = _tabulate(R0, identity, nothing)

# Durations of the steps of the construction, for the measurements of J1.c (c3): with a
# vector `timings`, each step pushes `name => seconds`; with `nothing`, nothing is done.
@inline _mark!(::Nothing, ::Symbol, t::UInt64) = t
function _mark!(timings::Vector{Pair{Symbol, Float64}}, name::Symbol, t::UInt64)
    now = time_ns()
    push!(timings, name => (now - t) / 1e9)
    return now
end

function _tabulate(R0::AbstractMatrix{<:Integer}, identity::Integer, timings)
    t = time_ns()
    N, n = size(R0)
    n >= 1 || _tab_error("no simple")
    n <= typemax(Int32) || _tab_error("$n simples, at most typemax(Int32) allowed")
    1 <= identity <= n || _tab_error("identity $identity out of 1:$n")
    e = Int32(identity)
    R = Matrix{Int32}(undef, N, n)
    for p in 1:n, a in 1:N
        v = R0[a, p]
        0 <= v <= n || _tab_error("entry R[$a, $p] = $v out of 0:$n")
        R[a, p] = v
    end
    atom = Vector{Int32}(undef, N)
    seen = falses(n)
    for a in 1:N
        x = R[a, e]
        x != 0 || _tab_error("atom $a missing: R[$a, identity] = 0")
        x != e || _tab_error("atom $a is the identity")
        seen[x] && _tab_error("atom $a repeated: simple $x")
        seen[x] = true
        atom[a] = x
    end
    t = _mark!(timings, :input, t)

    # 1. Lengths, by a breadth-first search from the identity; homogeneity.
    len = fill(Int32(-1), n)
    len[e] = 0
    bylen = Vector{Int32}(undef, n)               # simples by nondecreasing length
    bylen[1] = e
    head, tail = 1, 1
    while head <= tail
        x = bylen[head]
        head += 1
        for a in 1:N
            y = R[a, x]
            y == 0 && continue
            if len[y] < 0
                len[y] = len[x] + 1
                tail += 1
                bylen[tail] = y
            elseif len[y] != len[x] + 1
                _tab_error("not homogeneous: the cover $x → $y by atom $a goes from length " *
                           "$(len[x]) to length $(len[y])")
            end
        end
    end
    tail == n || _tab_error("simple $(findfirst(<(0), len)) not reached from the identity")
    maxlen = len[bylen[n]]
    maxlen <= 255 || _tab_error("length $maxlen of simple $(bylen[n]) larger than 255")
    t = _mark!(timings, :lengths, t)

    # Right cancellation: ua = u'a implies u = u', so the rows of R are injective.
    mark = zeros(Int32, n)
    for a in 1:N, x in 1:n
        y = R[a, x]
        y == 0 && continue
        mark[y] == a && _tab_error("right cancellation fails: two simples times atom $a give $y")
        mark[y] = a
    end
    t = _mark!(timings, :right_cancellation, t)

    # 2. Left covers by the atoms, L[a, u] = au, propagated along the right covers:
    # a(ub) = (au)b; the values reached by different paths must agree (associativity).
    L = zeros(Int32, N, n)
    done = falses(n)
    for a in 1:N
        L[a, e] = atom[a]
    end
    done[e] = true
    for u in bylen, b in 1:N
        v = R[b, u]
        v == 0 && continue
        first_path = !done[v]
        for a in 1:N
            au = L[a, u]
            w = au == 0 ? Int32(0) : R[b, au]
            if first_path
                L[a, v] = w
            elseif L[a, v] != w
                _tab_error("not associative: atom $a times simple $v (= $u · atom $b) " *
                           "gives $(L[a, v]) and $w along two paths")
            end
        end
        done[v] = true
    end
    t = _mark!(timings, :left_covers, t)

    # 3. Left quotients by the atoms (left cancellation), descents, boundedness.
    Q = zeros(Int32, N, n)
    for u in 1:n, a in 1:N
        ν = L[a, u]
        ν == 0 && continue
        Q[a, ν] == 0 || _tab_error("left cancellation fails: atom $a times $(Q[a, ν]) and " *
                                   "times $u give $ν")
        Q[a, ν] = u
    end
    t = _mark!(timings, :left_quotients, t)
    maximal = [x for x in 1:n if all(iszero, view(R, :, x))]
    length(maximal) == 1 ||
        _tab_error("not bounded: $(length(maximal)) simples without right cover ($(first(maximal, 5)))")
    firstdesc = zeros(Int32, n)                   # smallest left descent, 0 for e
    for ν in 1:n
        a = findfirst(!iszero, view(Q, :, ν))
        firstdesc[ν] = a === nothing ? 0 : a
    end
    t = _mark!(timings, :bounded, t)

    # 4. Canonical numbering, stratum by stratum: key (first letter, rank of the rest).
    canon = Vector{Int32}(undef, n)               # provisional -> canonical
    srank = zeros(Int32, n)                       # rank in its stratum
    canon[e] = 0
    srank[e] = 1
    lo = 2
    while lo <= n
        k = len[bylen[lo]]
        hi = lo
        while hi < n && len[bylen[hi + 1]] == k
            hi += 1
        end
        layer = view(bylen, lo:hi)
        layer_keys = [(firstdesc[ν], srank[Q[firstdesc[ν], ν]]) for ν in layer]
        for (r, i) in enumerate(sortperm(layer_keys))
            ν = layer[i]
            srank[ν] = r
            canon[ν] = lo - 2 + r
        end
        lo = hi + 1
    end
    t = _mark!(timings, :canonical, t)

    # 5. Lexicographic ranks: the smallest words are closed under prefixes, so they form
    # a tree, whose preorder (children by last letter) is the lexicographic order.
    parent = zeros(Int32, n)
    lastletter = zeros(Int32, n)
    for ν in view(bylen, 2:n)
        a = firstdesc[ν]
        rest = Q[a, ν]
        if rest == e
            parent[ν], lastletter[ν] = e, a
        else
            parent[ν] = L[a, parent[rest]]
            lastletter[ν] = lastletter[rest]
            parent[ν] != 0 || _tab_error("internal error: prefix of $ν not simple")
        end
    end
    children = sort!([(parent[ν], lastletter[ν], ν) for ν in 1:n if ν != e])
    start = zeros(Int, n + 1)                     # children of x: start[x]:start[x + 1] - 1
    for (p, _, _) in children
        start[p + 1] += 1
    end
    start[1] = 1
    for x in 1:n
        start[x + 1] += start[x]
    end
    lexrank = Vector{Int32}(undef, n)
    stack = Int32[e]
    r = Int32(0)
    while !isempty(stack)
        x = pop!(stack)
        lexrank[x] = r
        r += Int32(1)
        for k in (start[x + 1] - 1):-1:start[x]
            push!(stack, children[k][3])
        end
    end
    t = _mark!(timings, :lex_ranks, t)

    K = max(1, cld(N, 64))
    G = _tabulate_canonical(Val(K), N, n, R, Q, len, canon, lexrank, timings, t)
    return G, canon
end

# The tables in the canonical numbering, and step 6: the lcms with the atoms, by a
# search of the multiples of each simple by layers of increasing length.
function _tabulate_canonical(::Val{K}, N, n, R, Q, len, canon, lexrank, timings, t) where {K}
    lengths = Vector{UInt8}(undef, n)
    ranks = Vector{Int32}(undef, n)
    Rc = Matrix{Int32}(undef, N, n)
    Qc = Matrix{Int32}(undef, N, n)
    for x in 1:n
        c = canon[x] + 1
        lengths[c] = len[x]
        ranks[c] = lexrank[x]
        for a in 1:N
            Rc[a, c] = R[a, x] == 0 ? Int32(-1) : canon[R[a, x]]
            Qc[a, c] = Q[a, x] == 0 ? Int32(-1) : canon[Q[a, x]]
        end
    end
    descents = Vector{StaticBitSet{K}}(undef, n)
    for c in 1:n
        descents[c] = StaticBitSet{K}(a for a in 1:N if Qc[a, c] >= 0)
    end
    t = _mark!(timings, :tables, t)

    lcms = Matrix{NTuple{2, Int32}}(undef, N, n)
    visited = zeros(Int32, n)                     # visited[m + 1] == ω + 1: m met for ω
    served = zeros(Int, N)                        # layer where atom a was served, 0 if not
    layer = NTuple{2, Int32}[]                    # pairs (t, ωt)
    next = NTuple{2, Int32}[]
    for ω in Int32(0):Int32(n - 1)
        fill!(served, 0)
        nserved = 0
        for a in descents[ω + 1]                  # a ≼ ω: (ω, 0)
            lcms[a, ω + 1] = (ω, Int32(0))
            served[a] = -1
            nserved += 1
        end
        empty!(layer)
        push!(layer, (Int32(0), ω))
        visited[ω + 1] = ω + 1
        depth = 0
        while nserved < N
            depth += 1
            empty!(next)
            for (t, m) in layer, b in 1:N
                m2 = Rc[b, m + 1]
                (m2 < 0 || visited[m2 + 1] == ω + 1) && continue
                visited[m2 + 1] = ω + 1
                t2 = Rc[b, t + 1]
                t2 >= 0 || _tab_error("not a Garside germ: a right divisor of $(m2) is not simple")
                push!(next, (t2, m2))
            end
            isempty(next) && _tab_error("no common multiple of simple $ω and an atom")
            for (t2, m2) in next, a in descents[m2 + 1]
                if served[a] == 0
                    lcms[a, ω + 1] = (m2, t2)
                    served[a] = depth
                    nserved += 1
                elseif served[a] == depth
                    _tab_error("no lcm: simple $ω and atom $a have two minimal common " *
                               "multiples of length $(lengths[m2 + 1]), $(lcms[a, ω + 1][1]) and $m2")
                end
            end
            layer, next = next, layer
        end
    end
    _mark!(timings, :lcms, t)
    return TableGerm{K}(N, lengths, descents, Rc, Qc, lcms, ranks)
end

"""
    TableGerm(R, identity)
    TableGerm(G::AbstractGerm)

The tabulated germ given by its right covers (see [`tabulate`](@ref)), or the
tabulation of a germ `G` whose simples are all numbered (see `tabulate(G)`).
"""
TableGerm(R::AbstractMatrix{<:Integer}, identity::Integer) = first(tabulate(R, identity))
TableGerm(G::AbstractGerm) = first(tabulate(G))

# Right covers of a germ whose simples are all numbered, in the provisional numbering
# s + 1; `left = true` gives the left covers a·s, which are the right covers of Mᵒᵖ.
function _covers(G::AbstractGerm; left::Bool = false)
    S = simples(G)
    N = natoms(G)
    R = zeros(Int32, N, length(S))
    for s in S, a in 1:N
        v = left ? product(G, a, s) : product(G, s, a)
        R[a, s + 1] = v === nothing ? Int32(0) : v + Int32(1)
    end
    return R
end

"""
    tabulate(G::AbstractGerm) -> (T, renumber)

The [`TableGerm`](@ref) `T` of a germ `G` whose simples are all numbered (`simples(G)`
is defined), built from its right covers, and the numbering: `renumber[s + 1]` is the
number in `T` of the simple `s` of `G` (the identity if `G` is in the canonical order).
"""
tabulate(G::AbstractGerm) = tabulate(_covers(G), 1)

"""
    opposite(G::AbstractGerm)

Generic method for a germ whose simples are all numbered: the [`TableGerm`](@ref) of the
opposite monoid, built from the left covers `a·s` of `G`, which are the right covers of
`Mᵒᵖ`. It is recomputed at each call.
"""
opposite(G::AbstractGerm) = TableGerm(_covers(G; left = true), 1)

# The tables of a TableGerm are never modified after its construction.
_thread_safe(::TableGerm) = true

# Braid monoids with more simples than this are not tabulated by `germ` but given a lazy
# CoxeterGerm (J1.c, c4: tabulation in a few seconds up to B⁺(D₇) and B⁺(A₈), 12.9 s for
# B⁺(B₇), 154 s for B⁺(E₇), measured by j1/c4_threshold.jl of CatLab).
const _TABULATION_THRESHOLD = 5 * 10^5
