# Germs: the interface between a Garside monoid and the construction of its resolution.

"""
    AbstractGerm <: Gauss.Complexes.AbstractMonoid

Supertype of the germs of Garside monoids, in the sense of *Foundations of Garside
Theory* (Dehornoy, Digne, Godelle, Krammer, Michel, ch. VI): the simples with their
partial product, from which the construction of the resolution reads everything it
needs. A germ is the descriptor of the monoid of a [`MonoidAlgebra`](@ref): it fixes the
numbering of the simples, which are the labels of the differentials.

# Numbering

Simples are numbered by `Int32`: `0` is the identity and `1:N` are the atoms, in the
order of the germ (for a germ built from Chevie, the order of `B.atoms`), so that the
atom `i` is the simple `i`. A tabulated germ numbers all its simples in the canonical
order, by length, then by smallest word (lexicographic on the atoms); its last simple is
then Δ. A lazy germ numbers its simples as it meets them, after the identity and the
atoms. The order of the atoms is part of the germ.

# Contract

- *Homogeneity*: atoms have length 1, every reduced word of a simple `s` has `ℓ(s)`
  letters, and `ℓ(st) = ℓ(s) + ℓ(t)` when `st` is simple.
- *Bounded Garside germ*: for a simple `ω` and an atom `a`, the right lcm `ω ∨ a` exists
  and is simple; if `s ≼ t` are simples, then `s \\ t` is simple; the simples are the
  divisors of Δ.
- *Side*: `≼` is left divisibility and lcms are right lcms (`ω ≼ ω ∨ a`), the convention
  of the code of Gauss; that of Dehornoy–Lafont and Garnier is obtained by
  [`opposite`](@ref).
- *Mutability*: a query may update the memo of a lazy germ, so that concurrent queries
  are only safe for a germ documented as immutable.

# Interface

Required for every germ: [`natoms`](@ref), `length(G, s)` (a method of `Base.length`,
`ℓ(s)`), [`product`](@ref), [`ldiv`](@ref), [`left_descents`](@ref),
[`lcm_atom`](@ref), [`opposite`](@ref). Required for a germ whose simples are all
numbered (tabulated germ), optional otherwise: [`nsimples`](@ref) (the number of
simples, a constant, also defined by a lazy germ that knows it without enumerating) and
[`simples`](@ref). Derived, with a generic method that a germ may replace: `one(G)`,
`isone(G, s)`, [`atoms`](@ref), [`is_left_divisor`](@ref), [`word`](@ref),
[`isless_lex`](@ref).

Arguments that are simples accept any integer; simples are returned as `Int32`, lengths
and counts as `Int`.
"""
abstract type AbstractGerm <: AbstractMonoid end

# ------------------------------------------------------------------------------------
# Required functions

"""
    natoms(G::AbstractGerm) -> Int

The number `N` of atoms of `G`; the atoms are the simples `1:N`.
"""
function natoms end

"""
    nsimples(G::AbstractGerm) -> Int

The number of simples of `G`, identity included: a constant. Defined by every tabulated
germ, whose simples are then `0:nsimples(G)-1`, and by a lazy germ only if it knows the
number without enumerating the simples.
"""
function nsimples end

"""
    simples(G::AbstractGerm)

The simples of `G`, `Int32(0):Int32(nsimples(G) - 1)`. Defined only by a germ whose
simples are all numbered (a tabulated germ), never by a lazy germ.
"""
function simples end

"""
    product(G::AbstractGerm, s, t) -> Union{Int32, Nothing}

The simple `st` if the product of the simples `s` and `t` is simple, `nothing`
otherwise. Where the product is known to exist, write `product(G, s, t)::Int32`: the
assertion documents the invariant and checks it.
"""
function product end

"""
    ldiv(G::AbstractGerm, s, t) -> Union{Int32, Nothing}

The simple `s \\ t`, that is the `u` such that `su = t`, if `s ≼ t`, `nothing`
otherwise.
"""
function ldiv end

"""
    left_descents(G::AbstractGerm, s) -> AbstractSet{Int}

The atoms `a` such that `a ≼ s`, iterated in increasing order, with a fast `in`. The
germs of Gauss return a [`StaticBitSet`](@ref).
"""
function left_descents end

"""
    lcm_atom(G::AbstractGerm, ω, a) -> NamedTuple{(:lcm, :complement), NTuple{2, Int32}}

The right lcm `ω ∨ a` of the simple `ω` and the atom `a`, and its complement `τ`, such
that `ωτ = ω ∨ a`: `(lcm = ω ∨ a, complement = τ)`, which destructures as
`ωa, τ = lcm_atom(G, ω, a)`. If `a ≼ ω`, the result is `(lcm = ω, complement = 0)`.
"""
function lcm_atom end

"""
    opposite(G::AbstractGerm)

The germ of the opposite monoid `Mᵒᵖ`: the same atoms, in the same order, and the
reversed product, `s ·ᵒᵖ t = ts`. Its left divisibility is the right divisibility of
`M`, its lcms are the left lcms of `M`. A tabulated opposite numbers its simples in its
own canonical order, so that `opposite(opposite(G))` numbers them as `G`.
"""
function opposite end

# ------------------------------------------------------------------------------------
# Derived functions

"""
    one(G::AbstractGerm) -> Int32

The identity of `G`, the simple `0`.
"""
Base.one(::AbstractGerm) = Int32(0)

"""
    isone(G::AbstractGerm, s) -> Bool

Whether the simple `s` is the identity.
"""
Base.isone(::AbstractGerm, s::Integer) = iszero(s)

"""
    atoms(G::AbstractGerm)

The atoms of `G`, `Int32(1):Int32(natoms(G))`.
"""
atoms(G::AbstractGerm) = Int32(1):Int32(natoms(G))

"""
    is_left_divisor(G::AbstractGerm, s, t) -> Bool

Whether `s ≼ t`. The generic method answers without division when `s` is the identity,
equal to `t`, an atom (by [`left_descents`](@ref)) or not shorter than `t`, and
otherwise asks [`ldiv`](@ref).
"""
function is_left_divisor(G::AbstractGerm, s::Integer, t::Integer)
    iszero(s) && return true
    s == t && return true
    s <= natoms(G) && return Int(s) in left_descents(G, t)
    length(G, s) < length(G, t) || return false
    return ldiv(G, s, t) !== nothing
end

"""
    word(G::AbstractGerm, s) -> Vector{Int32}

The smallest reduced word of the simple `s` in the atoms, for the lexicographic order:
at each step, the smallest left descent. Empty for the identity.
"""
function word(G::AbstractGerm, s::Integer)
    n = length(G, s)
    w = Vector{Int32}(undef, n)
    r = Int32(s)
    for k in 1:n
        D = left_descents(G, r)
        isempty(D) && throw(ArgumentError("simple $r of length > 0 without left descent"))
        a = Int32(first(D))
        w[k] = a
        r = ldiv(G, a, r)::Int32
    end
    iszero(r) || throw(ArgumentError("length($s) = $n, but $n left divisions leave the simple $r"))
    return w
end

"""
    isless_lex(G::AbstractGerm, s, t) -> Bool

Whether the smallest word of `s` comes strictly before that of `t` in the lexicographic
order, a prefix coming first. The generic method compares the words letter by letter,
without allocating; a tabulated germ compares precomputed ranks.
"""
function isless_lex(G::AbstractGerm, s::Integer, t::Integer)
    x, y = Int32(s), Int32(t)
    while true
        x == y && return false
        iszero(y) && return false
        iszero(x) && return true
        a = Int32(first(left_descents(G, x)))
        b = Int32(first(left_descents(G, y)))
        a == b || return a < b
        x = ldiv(G, a, x)::Int32
        y = ldiv(G, b, y)::Int32
    end
end

# ------------------------------------------------------------------------------------
# Entry point

"""
    germ(M)

The germ of the monoid `M`, of the type best suited to its size. Methods are added by
the package extensions: with Chevie loaded, `M` may be a `BraidMonoid` of a finite
Coxeter group or a `DualBraidMonoid` (of a Coxeter group, or of a well-generated
irreducible complex reflection group), and the result is a [`TableGerm`](@ref) if the
number of simples, computed from the degrees, is at most 10⁵. Above, an `ArgumentError`
is thrown for now: the lazy germ of the large classical monoids is still to be written;
`TableGerm(M)` forces the tabulation.
"""
function germ end
