"""
    GaussChevieExt

Package extension loaded together with Chevie: the germs of the braid monoids
(`BraidMonoid` of a finite Coxeter group) and of the dual braid monoids
(`DualBraidMonoid`) of Chevie, built by `Gauss.germ(M)` and `Gauss.TableGerm(M)`.

The right covers σ → σa given to `Gauss.tabulate` come from the arithmetic of the group
alone, without lcms of Chevie: for a braid monoid, by a search of the Cayley graph of the
generators; for a dual braid monoid, from the simples listed by `left_divisors`. An
element of the group is identified by the images of the roots of the generators, which
form a basis: a key of `n` small integers, updated by the action of the atoms without
forming the product of permutations.
"""
module GaussChevieExt

import Gauss
import Chevie

const ChevieMonoid = Union{Chevie.BraidMonoid, Chevie.DualBraidMonoid}

"""
    _simples_count(M) -> BigInt

The number of simples of `M`, from the degrees of its group, without enumeration: |W|
for a braid monoid; for a dual braid monoid of an irreducible group, the Catalan number
Π (dᵢ + h) / dᵢ, `h` being the largest degree, and |W| (an upper bound) otherwise.
"""
_simples_count(M::Chevie.BraidMonoid) = prod(BigInt, Chevie.degrees(M.W))

function _simples_count(M::Chevie.DualBraidMonoid)
    d = Chevie.degrees(M.W)
    length(Chevie.refltype(M.W)) == 1 || return prod(BigInt, d)
    h = maximum(d)
    return prod(BigInt, d .+ h) ÷ prod(BigInt, d)
end

function Gauss.germ(M::ChevieMonoid)
    n = _simples_count(M)
    n <= Gauss._TABULATION_THRESHOLD ||
        throw(ArgumentError("germ: $M has $n simples, above the tabulation threshold " *
                            "$(Gauss._TABULATION_THRESHOLD); the lazy germ of the large " *
                            "classical monoids is not written yet; TableGerm(M) forces " *
                            "the tabulation"))
    return Gauss.TableGerm(M)
end

Gauss.TableGerm(M::ChevieMonoid) = Gauss.TableGerm(_right_covers(M)...)

# ------------------------------------------------------------------------------------
# Right covers

# Images of the points 1:D under each atom, as vectors.
_atom_images(atoms, D) = [Int16[Int16(i^a) for i in 1:D] for a in atoms]

# The key of an element w: the images of the roots of the generators, 1:n.
_key(w, ::Val{n}) where {n} = ntuple(j -> Int16(j^w), Val(n))

# The key of σa from that of σ, by i^(σa) = (i^σ)^a.
_shift(k::NTuple{n, Int16}, image::Vector{Int16}) where {n} = map(x -> image[x], k)

_perm_degree(p) = length(Chevie.perm(p))

"""
    _right_covers(M, timings = nothing) -> (R, identity)

The right covers of the simples of `M`, in the form taken by `Gauss.tabulate`, the
identity being numbered 1. With a vector `timings`, the durations of the enumeration of
the simples (dual monoids only) and of the covers are pushed, for the measurements of
J1.c.
"""
_right_covers(M::ChevieMonoid, timings = nothing) =
    _right_covers(M, Val(length(Chevie.gens(M.W))), timings)

# Braid monoid: breadth-first search of the Cayley graph of the generators. For σ of
# length k and a generator s, σs has length k ± 1: it is a cover if it is new or of
# length k + 1, and s is a right descent of σ otherwise.
function _right_covers(M::Chevie.BraidMonoid, ::Val{nb}, timings) where {nb}
    t = time_ns()
    atoms = M.atoms
    N = length(atoms)
    total = _simples_count(M)
    total <= typemax(Int32) || throw(ArgumentError("TableGerm: $M has $total simples, too many"))
    n = Int(total)
    D = max(nb, maximum(_perm_degree, atoms; init = 0))
    images = _atom_images(atoms, D)
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

# Dual braid monoid: the simples by length from left_divisors, then σa is a cover if it
# is a simple of length ℓ(σ) + 1.
function _right_covers(M::Chevie.DualBraidMonoid, ::Val{nb}, timings) where {nb}
    t = time_ns()
    layers = Chevie.left_divisors(M, M.δ)
    t = Gauss._mark!(timings, :enumeration, t)
    simples = reduce(vcat, layers)
    n = length(simples)
    n <= typemax(Int32) || throw(ArgumentError("TableGerm: $M has $n simples, too many"))
    isone(simples[1]) || throw(ArgumentError("TableGerm: the first simple of $M is not the identity"))
    lengths = Int32[k - 1 for (k, l) in enumerate(layers) for _ in l]
    atoms = M.atoms
    N = length(atoms)
    D = max(nb, maximum(_perm_degree, atoms; init = 0), maximum(_perm_degree, simples))
    images = _atom_images(atoms, D)
    simple_keys = [_key(w, Val(nb)) for w in simples]
    index = Dict{NTuple{nb, Int16}, Int32}(k => Int32(p) for (p, k) in enumerate(simple_keys))
    length(index) == n ||
        throw(ArgumentError("TableGerm: $(length(index)) distinct keys for the $n simples of " *
                            "$M: the images of the roots of the generators do not determine " *
                            "the elements"))
    R = zeros(Int32, N, n)
    for σ in 1:n, a in 1:N
        τ = get(index, _shift(simple_keys[σ], images[a]), Int32(0))
        τ != 0 && lengths[τ] == lengths[σ] + 1 && (R[a, σ] = τ)
    end
    Gauss._mark!(timings, :covers, t)
    return R, 1
end

end # module GaussChevieExt
