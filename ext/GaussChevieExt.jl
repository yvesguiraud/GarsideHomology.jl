"""
    GaussChevieExt

Package extension loaded together with Chevie: the germs of the braid monoids
(`BraidMonoid` of a finite Coxeter group) and of the dual braid monoids of Chevie,
built by `Gauss.germ(M)` and `Gauss.TableGerm(M)`, and `Gauss.dual_germ(W; c)`, which
builds the germ of a dual braid monoid from the group without building the monoid.

The right covers σ → σa given to `Gauss.tabulate` come from the arithmetic of the group
alone. An element is identified by the images of the roots of the generators, which
form a basis: a key of `n` small integers, updated by the action of the atoms without
forming products of permutations.
- Braid monoid: search of the Cayley graph of the generators.
- Dual braid monoid of a real group (2 among the degrees): search by layers with the
  criterion of Brady and Watt, σa is simple if and only if the root of a lies in the
  moved space of the Kreweras complement σ⁻¹c, by exact linear algebra.
- Dual braid monoid of a complex group: prefixes of the Hurwitz orbit of a reduced
  decomposition of c (Bessis), then covers by the layers.
- Reference: the simples listed by `left_divisors` of Chevie (`method = :chevie`).
"""
module GaussChevieExt

import Gauss
import Chevie

const ChevieMonoid = Union{Chevie.BraidMonoid, Chevie.DualBraidMonoid}
const ReflectionGroup = Union{Chevie.FiniteCoxeterGroup, Chevie.PermRootGroup}

# ------------------------------------------------------------------------------------
# Numbers of simples, threshold

_order(W) = prod(BigInt, Chevie.degrees(W))

"""
    _simples_count(M) -> BigInt

The number of simples of `M`, from the degrees of its group, without enumeration: |W|
for a braid monoid; for a dual braid monoid of an irreducible group, the Catalan number
Π (dᵢ + h) / dᵢ, `h` being the largest degree, and |W| (an upper bound) otherwise.
"""
_simples_count(M::Chevie.BraidMonoid) = _order(M.W)
_simples_count(M::Chevie.DualBraidMonoid) = _dual_count(M.W)

function _dual_count(W)
    d = Chevie.degrees(W)
    length(Chevie.refltype(W)) == 1 || return prod(BigInt, d)
    h = maximum(d)
    return prod(BigInt, d .+ h) ÷ prod(BigInt, d)
end

function _check_threshold(n, what)
    n <= Gauss._TABULATION_THRESHOLD ||
        throw(ArgumentError("$what has $n simples, above the tabulation threshold " *
                            "$(Gauss._TABULATION_THRESHOLD), and there is no lazy germ " *
                            "for the dual braid monoids; TableGerm forces the tabulation"))
    return nothing
end

function Gauss.germ(M::Chevie.DualBraidMonoid)
    _check_threshold(_simples_count(M), "germ: $M")
    return Gauss.TableGerm(M)
end

# Braid monoid: tabulated under the threshold, lazy above.
Gauss.germ(M::Chevie.BraidMonoid) =
    _simples_count(M) <= Gauss._TABULATION_THRESHOLD ? Gauss.TableGerm(M) :
    Gauss.CoxeterGerm(Chevie.cartan(M.W))

Gauss.TableGerm(M::ChevieMonoid) = Gauss.TableGerm(_right_covers(M)...)

function Gauss.dual_germ(W::ReflectionGroup; c = nothing)
    n = length(Chevie.gens(W))
    length(Chevie.refltype(W)) == 1 && n == Chevie.semisimplerank(W) ||
        throw(ArgumentError("dual_germ: $W must be well-generated and irreducible"))
    _check_threshold(_dual_count(W), "dual_germ: $W")
    word = c === nothing ? _default_coxeter_word(W, Val(n)) : collect(Int, c)
    return Gauss.TableGerm(_dual_covers(W, W(word...), word, nothing)...)
end

# ------------------------------------------------------------------------------------
# Keys

# Images of the points 1:D under each element, as vectors.
_images(elements, D) = [Int16[Int16(i^a) for i in 1:D] for a in elements]

# The key of an element w: the images of the roots of the generators, 1:n.
_key(w, ::Val{n}) where {n} = ntuple(j -> Int16(j^w), Val(n))

# The key of σa from that of σ, by i^(σa) = (i^σ)^a.
_shift(k::NTuple{n, Int16}, image::Vector{Int16}) where {n} = map(x -> image[x], k)

_perm_degree(p) = length(Chevie.perm(p))

# ------------------------------------------------------------------------------------
# Exact linear algebra over a field (rationals or cyclotomics), on small matrices

# A basis of the vectors x with B x = 0.
function _nullspace(B::AbstractMatrix)
    T = eltype(B)
    M = Matrix{T}(B)
    m, n = size(M)
    pivots = Int[]
    r = 0
    for col in 1:n
        r == m && break
        p = findfirst(i -> !iszero(M[i, col]), (r + 1):m)
        p === nothing && continue
        p += r
        r += 1
        if p != r
            for j in 1:n
                M[p, j], M[r, j] = M[r, j], M[p, j]
            end
        end
        piv = M[r, col]
        for j in 1:n
            M[r, j] = M[r, j] / piv
        end
        for i in 1:m
            i == r && continue
            f = M[i, col]
            iszero(f) && continue
            for j in 1:n
                M[i, j] -= f * M[r, j]
            end
        end
        push!(pivots, col)
    end
    basis = Vector{T}[]
    for free in setdiff(1:n, pivots)
        x = [zero(T) for _ in 1:n]
        x[free] = one(T)
        for (i, pc) in enumerate(pivots)
            x[pc] = -M[i, free]
        end
        push!(basis, x)
    end
    return basis
end

_rank(A::AbstractMatrix) = size(A, 2) - length(_nullspace(A))

_minus_identity(X) = [X[i, j] - (i == j ? one(eltype(X)) : zero(eltype(X))) for i in axes(X, 1), j in axes(X, 2)]

"""
    _Roots

The coordinates of the roots of a reflection group in the basis of the roots of its
generators (numbers 1:n), in a field `T` (rationals, or cyclotomics), and as integers
when they are all integral (crystallographic case).
"""
struct _Roots{T}
    coords::Vector{Vector{T}}
    integral::Bool
end

_rationalize(x) = Chevie.conductor(x) == 1 ? Rational{Int}(x) : x

function _Roots(W, n)
    roots = Chevie.roots(W)
    base = permutedims(reduce(hcat, roots[1:n]))
    raw = [Chevie.solutionmat(base, v) for v in roots]
    any(isnothing, raw) && throw(ArgumentError("the roots of $W are not in the span of those of its generators"))
    values = map(v -> map(_rationalize, v), raw)
    if all(v -> all(x -> x isa Rational, v), values)
        coords = [Rational{Int}.(v) for v in values]
        return _Roots{Rational{Int}}(coords, all(v -> all(isinteger, v), coords))
    end
    T = promote_type(map(eltype, raw)...)
    return _Roots{T}([T.(v) for v in raw], false)
end

# The matrix of w in the basis of the roots of the generators, from its key: column j is
# the coordinate vector of the root j^w.
_matrix(R::_Roots, k::NTuple) = reduce(hcat, [R.coords[i] for i in k])

# A vector of integers proportional to a vector of rationals.
function _integral(v::Vector{Rational{Int}})
    d = foldl(lcm, (denominator(x) for x in v); init = 1)
    return [Int(x * d) for x in v]
end

# ------------------------------------------------------------------------------------
# Covers from the simples (keys by length): σa is a cover if it is a simple of length
# ℓ(σ) + 1.

function _covers_by_layers(simple_keys::Vector{NTuple{nb, Int16}}, lengths, images) where {nb}
    n = length(simple_keys)
    index = Dict{NTuple{nb, Int16}, Int32}(k => Int32(p) for (p, k) in enumerate(simple_keys))
    length(index) == n ||
        throw(ArgumentError("TableGerm: $(length(index)) distinct keys for $n simples: the " *
                            "images of the roots of the generators do not determine the elements"))
    N = length(images)
    R = zeros(Int32, N, n)
    for σ in 1:n, a in 1:N
        τ = get(index, _shift(simple_keys[σ], images[a]), Int32(0))
        τ != 0 && lengths[τ] == lengths[σ] + 1 && (R[a, σ] = τ)
    end
    return R
end

# ------------------------------------------------------------------------------------
# Right covers

"""
    _right_covers(M, timings = nothing; method = :auto) -> (R, identity)

The right covers of the simples of `M`, in the form taken by `Gauss.tabulate`, the
identity being numbered 1. For a dual braid monoid, `method = :chevie` uses the simples
listed by `left_divisors` (the reference of J1.c, c3); `:auto` uses the criterion of
Brady and Watt for a real group and the Hurwitz orbit otherwise. With a vector
`timings`, the durations of the enumeration of the simples and of the covers are pushed.
"""
_right_covers(M::Chevie.BraidMonoid, timings = nothing) =
    _classical_covers(M, Val(length(Chevie.gens(M.W))), timings)

function _right_covers(M::Chevie.DualBraidMonoid, timings = nothing; method::Symbol = :auto)
    method === :chevie && return _chevie_dual_covers(M, Val(length(Chevie.gens(M.W))), timings)
    method === :auto || throw(ArgumentError("unknown method $method"))
    return _dual_covers(M.W, M.δ, M.δword, timings; atoms = M.atoms)
end

# Braid monoid: breadth-first search of the Cayley graph of the generators. For σ of
# length k and a generator s, σs has length k ± 1: it is a cover if it is new or of
# length k + 1, and s is a right descent of σ otherwise.
function _classical_covers(M::Chevie.BraidMonoid, ::Val{nb}, timings) where {nb}
    t = time_ns()
    atoms = M.atoms
    N = length(atoms)
    total = _simples_count(M)
    total <= typemax(Int32) || throw(ArgumentError("TableGerm: $M has $total simples, too many"))
    n = Int(total)
    D = max(nb, maximum(_perm_degree, atoms; init = 0))
    images = _images(atoms, D)
    R = zeros(Int32, N, n)
    lengths = Vector{Int32}(undef, n)
    simple_keys = Vector{NTuple{nb, Int16}}(undef, n)
    index = Dict{NTuple{nb, Int16}, Int32}()
    sizehint!(index, n)
    simple_keys[1] = ntuple(j -> Int16(j), Val(nb))
    index[simple_keys[1]] = 1
    lengths[1] = 0
    found = 1
    head = 1
    while head <= found
        kσ = simple_keys[head]
        for a in 1:N
            kτ = _shift(kσ, images[a])
            τ = get(index, kτ, Int32(0))
            if τ == 0
                found < n || throw(ArgumentError("TableGerm: more than $n elements in $(M.W): " *
                                                 "the degrees are wrong"))
                found += 1
                τ = Int32(found)
                simple_keys[τ] = kτ
                index[kτ] = τ
                lengths[τ] = lengths[head] + 1
                R[a, head] = τ
            elseif lengths[τ] == lengths[head] + 1
                R[a, head] = τ
            end
        end
        head += 1
    end
    found == n || throw(ArgumentError("TableGerm: $found distinct keys for the $n elements of " *
                                      "$(M.W): the images of the roots of the generators do " *
                                      "not determine the elements"))
    Gauss._mark!(timings, :covers, t)
    return R, 1
end

# Reference for a dual braid monoid: the simples by length from left_divisors.
function _chevie_dual_covers(M::Chevie.DualBraidMonoid, ::Val{nb}, timings) where {nb}
    t = time_ns()
    layers = Chevie.left_divisors(M, M.δ)
    t = Gauss._mark!(timings, :enumeration, t)
    simples = reduce(vcat, layers)
    length(simples) <= typemax(Int32) || throw(ArgumentError("TableGerm: $M has too many simples"))
    isone(simples[1]) || throw(ArgumentError("TableGerm: the first simple of $M is not the identity"))
    lengths = Int32[k - 1 for (k, l) in enumerate(layers) for _ in l]
    D = max(nb, maximum(_perm_degree, M.atoms; init = 0), maximum(_perm_degree, simples))
    R = _covers_by_layers([_key(w, Val(nb)) for w in simples], lengths, _images(M.atoms, D))
    Gauss._mark!(timings, :covers, t)
    return R, 1
end

# Dual braid monoid from the group W, the Coxeter element c (a permutation) and a word of
# c in the generators. `atoms`: those of a monoid of Chevie, or nothing (then the
# reflections ≤ c, in the order of unique(refls(W)), as in Chevie).
function _dual_covers(W, c, word, timings; atoms = nothing)
    nb = length(Chevie.gens(W))
    return _dual_covers(W, c, word, Val(nb), timings, atoms)
end

function _dual_covers(W, c, word, ::Val{nb}, timings, atoms) where {nb}
    t = time_ns()
    roots = _Roots(W, nb)
    reflections = unique(Chevie.refls(W))
    expected = _dual_count(W)
    if 2 in Chevie.degrees(W)
        atomlist = atoms === nothing ? reflections : atoms
        length(atomlist) == length(reflections) ||
            throw(ArgumentError("TableGerm: $(length(atomlist)) atoms for a real group with " *
                                "$(length(reflections)) reflections"))
        R, n = _real_covers(roots, c, atomlist, Val(nb))
        t = Gauss._mark!(timings, :enumeration, t)
        Gauss._mark!(timings, :covers, t)
    else
        simple_keys, lengths, atomlist = _hurwitz_simples(W, roots, c, word, reflections, atoms, Val(nb))
        t = Gauss._mark!(timings, :enumeration, t)
        D = max(nb, maximum(_perm_degree, atomlist; init = 0))
        R = _covers_by_layers(simple_keys, lengths, _images(atomlist, D))
        n = length(simple_keys)
        Gauss._mark!(timings, :covers, t)
    end
    length(Chevie.refltype(W)) == 1 && n != expected &&
        throw(ArgumentError("TableGerm: $n simples for $W, Π (dᵢ + h)/dᵢ = $expected"))
    return R, 1
end

# ------------------------------------------------------------------------------------
# Real groups: search by layers with the criterion of Brady and Watt. Each simple σ
# carries the matrix of its Kreweras complement σ̄ = σ⁻¹c; σa is a cover if and only if
# the root of a lies in Mov(σ̄) = Im(σ̄ − 1), and then the complement of σa is aσ̄.

function _real_covers(roots::_Roots{T}, c, atoms, ::Val{nb}) where {T, nb}
    N = length(atoms)
    D = max(nb, maximum(_perm_degree, atoms; init = 0), _perm_degree(c))
    images = _images(atoms, D)
    Xa = [_matrix(roots, _key(a, Val(nb))) for a in atoms]
    α = map(Xa) do X
        Y = _minus_identity(X)
        j = findfirst(j -> !all(iszero, view(Y, :, j)), 1:nb)
        j === nothing && throw(ArgumentError("TableGerm: an atom acts trivially"))
        Y[:, j]
    end
    if roots.integral
        return _real_search(Int.(_matrix(roots, _key(c, Val(nb)))), [Int.(X) for X in Xa],
                            [Int.(a) for a in α], images, Val(nb), true)
    else
        return _real_search(_matrix(roots, _key(c, Val(nb))), Xa, α, images, Val(nb), false)
    end
end

# The test vectors of Mov(X): a basis of the rows u with u (X − 1) = 0, as integers on
# the integral path.
function _test_vectors(X, integral::Bool)
    Y = _minus_identity(X)
    if integral
        return [_integral(u) for u in _nullspace(permutedims(Rational{Int}.(Y)))]
    else
        return _nullspace(permutedims(Y))
    end
end

function _real_search(Xc::Matrix{T}, Xa, α, images, ::Val{nb}, integral) where {T, nb}
    N = length(Xa)
    simple_keys = NTuple{nb, Int16}[ntuple(j -> Int16(j), Val(nb))]
    complement = Dict{Int32, Matrix{T}}(Int32(1) => Xc)
    index = Dict{NTuple{nb, Int16}, Int32}(simple_keys[1] => 1)
    covers = Tuple{Int, Int32, Int32}[]                  # (a, σ, σa)
    layer = Int32[1]
    k = 0
    while !isempty(layer)
        next = Int32[]
        for σ in layer
            X = pop!(complement, σ)
            U = _test_vectors(X, integral)
            length(U) == k || throw(ArgumentError("TableGerm: the complement of a simple of " *
                "length $k has a moved space of dimension $(nb - length(U)), not $(nb - k) " *
                "(lemma of Carter): the group is not real"))
            for a in 1:N
                all(u -> iszero(sum(u[i] * α[a][i] for i in 1:nb)), U) || continue
                kτ = _shift(simple_keys[σ], images[a])
                τ = get(index, kτ, Int32(0))
                if τ == 0
                    push!(simple_keys, kτ)
                    τ = Int32(length(simple_keys))
                    index[kτ] = τ
                    complement[τ] = X * Xa[a]
                    push!(next, τ)
                end
                push!(covers, (a, σ, τ))
            end
        end
        layer = next
        k += 1
    end
    n = length(simple_keys)
    R = zeros(Int32, N, n)
    for (a, σ, τ) in covers
        R[a, σ] = τ
    end
    return R, n
end

# ------------------------------------------------------------------------------------
# Complex groups: Hurwitz orbit of a reduced decomposition of c.

# A reduced decomposition of c into n reflections: the word in the generators if it has
# n letters, otherwise a search in depth pruned by the codimension of the fixed space
# (a product of n reflections equal to c is reduced, since ℓ_R(c) = n).
function _reduced_decomposition(W, roots, c, word, reflections, ::Val{nb}) where {nb}
    gens = Chevie.gens(W)
    if length(word) == nb
        d = [gens[i] for i in word]
        all(in(Set(reflections)), d) && prod(d) == c && return d
    end
    rset = Set(reflections)
    codim(w) = _rank(_minus_identity(_matrix(roots, _key(w, Val(nb)))))
    codim(c) == nb || throw(ArgumentError("TableGerm: c has a fixed space of codimension " *
                                          "$(codim(c)), not $nb"))
    function search(target, k)
        k == 1 && return target in rset ? [target] : nothing
        for r in reflections
            rest = inv(r) * target                      # target = r · rest
            codim(rest) == k - 1 || continue
            tail = search(rest, k - 1)
            tail === nothing || return vcat([r], tail)
        end
        return nothing
    end
    d = search(c, nb)
    d === nothing && throw(ArgumentError("TableGerm: no decomposition of c into $nb reflections"))
    return d
end

function _hurwitz_simples(W, roots, c, word, reflections, atoms, ::Val{nb}) where {nb}
    decomposition = _reduced_decomposition(W, roots, c, word, reflections, Val(nb))
    letters = eltype(reflections)[]
    number = Dict{eltype(reflections), Int32}()
    function letter(p)
        get!(number, p) do
            push!(letters, p)
            Int32(length(letters))
        end
    end
    start = ntuple(i -> letter(decomposition[i]), Val(nb))
    table = Dict{Tuple{Int32, Int32}, Int32}()
    conjugate(i, j) = get!(() -> letter(letters[i] * letters[j] * inv(letters[i])), table, (i, j))
    orbit = Gauss.hurwitz_orbit(start, conjugate)
    d = Chevie.degrees(W)
    expected = factorial(big(nb)) * big(maximum(d))^nb ÷ _order(W)
    length(orbit) == expected ||
        throw(ArgumentError("TableGerm: Hurwitz orbit of $(length(orbit)) decompositions, " *
                            "n!·hⁿ/|W| = $expected"))
    letterset = Set(letters)
    atomlist = atoms === nothing ? [r for r in reflections if r in letterset] : atoms
    length(atomlist) == length(letters) && Set(atomlist) == letterset ||
        throw(ArgumentError("TableGerm: the letters of the Hurwitz orbit are not the atoms"))
    D = max(nb, maximum(_perm_degree, letters; init = 0))
    images = _images(letters, D)
    lengthof = Dict{NTuple{nb, Int16}, Int32}()
    e = ntuple(j -> Int16(j), Val(nb))
    lengthof[e] = 0
    for t in orbit
        k = e
        for j in 1:nb
            k = _shift(k, images[t[j]])
            get!(lengthof, k, Int32(j))
        end
    end
    simple_keys = sort!(collect(keys(lengthof)); by = k -> (lengthof[k], k))
    return simple_keys, Int32[lengthof[k] for k in simple_keys], atomlist
end

# ------------------------------------------------------------------------------------
# Springer's criterion: c is a Coxeter element if its eigenspace for e^{2iπ/h} is a line
# (a(h) = 1) not contained in any reflecting hyperplane.

function _is_regular(W, roots, c, reflections, ::Val{nb}) where {nb}
    h = maximum(Chevie.degrees(W))
    ζ = Chevie.E(h)
    X = _matrix(roots, _key(c, Val(nb)))
    M = [X[i, j] - (i == j ? ζ : 0 * ζ) for i in 1:nb, j in 1:nb]
    V = _nullspace(M)
    length(V) == 1 || return false
    v = V[1]
    for r in reflections
        Xr = _matrix(roots, _key(r, Val(nb)))
        all(i -> iszero(sum(Xr[i, j] * v[j] for j in 1:nb) - v[i]), 1:nb) && return false
    end
    return true
end

# The next permutation of p in the lexicographic order, in place; false after the last.
function _next_permutation!(p::Vector{Int})
    i = findlast(k -> p[k] < p[k + 1], 1:(length(p) - 1))
    i === nothing && return false
    j = findlast(k -> p[k] > p[i], eachindex(p))
    p[i], p[j] = p[j], p[i]
    reverse!(p, i + 1, length(p))
    return true
end

# The default Coxeter element of dual_germ: the product 1:n of the generators if it is
# regular, otherwise the first order of the generators, in the lexicographic order, that
# is (with a message naming it); an error if there is none.
function _default_coxeter_word(W, ::Val{nb}) where {nb}
    roots = _Roots(W, nb)
    reflections = unique(Chevie.refls(W))
    p = collect(1:nb)
    while true
        if _is_regular(W, roots, W(p...), reflections, Val(nb))
            p == 1:nb || @info "dual_germ: 1:$nb is not a Coxeter element of $W; c = $p"
            return p
        end
        _next_permutation!(p) || break
    end
    throw(ArgumentError("dual_germ: no order of the generators of $W is a Coxeter element " *
                        "(Springer's criterion); give c"))
end

end # module GaussChevieExt
