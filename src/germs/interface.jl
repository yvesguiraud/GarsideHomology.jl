# Germs: the interface between a Garside monoid and the construction of its resolution.

"""
    AbstractGerm <: GarsideHomology.Complexes.AbstractMonoid

Supertype of the germs of Garside monoids, in the sense of *Foundations of Garside Theory* (Dehornoy, Digne, Godelle, Krammer, Michel, ch. VI). A germ is given by the simples and their partial product, from which the construction of the resolution reads everything it needs. A germ is also the descriptor of the monoid of a [`MonoidAlgebra`](@ref): it fixes the numbering of the simples, which are the labels of the differentials.

# Numbering

Simples are numbered by `Int32`: `0` is the identity, and `1:N` are the atoms, in the order of the germ, so that the atom `i` is the simple `i`. For a germ built from a monoid `M` of Chevie, the order of the atoms is that of `M.atoms`. The order of the atoms is part of the germ.
- A tabulated germ numbers all its simples in the canonical order: by length, then by smallest word, the words being compared in the lexicographic order on the atoms. Its last simple is then Δ.
- A lazy germ numbers its simples as it meets them, after the identity and the atoms.

# Contract

- *Homogeneity*: the atoms have length 1, every reduced word of a simple `s` has `ℓ(s)` letters, and `ℓ(st) = ℓ(s) + ℓ(t)` when `st` is simple.
- *Bounded Garside germ*: for a simple `ω` and an atom `a`, the right lcm `ω ∨ a` exists and is simple; if `s ≼ t` are simples, then `s \\ t` is simple; the simples are the divisors of Δ.
- *Side*: `≼` is left divisibility and the lcms are right lcms (`ω ≼ ω ∨ a`), which is the convention of GarsideHomology; the convention of Dehornoy–Lafont and Garnier is obtained by [`opposite`](@ref).
- *Mutability*: a query may update the memo of a lazy germ, so that concurrent queries are safe only for a germ documented as immutable.

# Interface

The functions of the interface fall into three groups:
- required for every germ: [`natoms`](@ref), `length(G, s)` (a method of `Base.length`, which gives `ℓ(s)`), [`product`](@ref), [`ldiv`](@ref), [`left_descents`](@ref), [`lcm_atom`](@ref) and [`opposite`](@ref);
- required for a germ whose simples are all numbered (a tabulated germ), and optional otherwise: [`nsimples`](@ref), and the internal method `GarsideHomology._simple_numbers(G)`, which returns the range `Int32(0):Int32(nsimples(G) - 1)` and from which [`simples`](@ref) and [`elements`](@ref) are derived;
- derived, with a generic method that a germ may replace: `one(G)`, `isone(G, s)`, [`atoms`](@ref), [`is_left_divisor`](@ref), [`word`](@ref), [`isless_lex`](@ref) and [`isless_deglex`](@ref).

The arguments that are simples accept any integer. Simples are returned as `Int32`, and lengths and counts as `Int`. For the user, a germ is also callable, and gives its elements as [`GermElement`](@ref)s, written as words in the atoms: `G(1, 2, 1)`, `G("121")`, `G()`.
"""
abstract type AbstractGerm <: AbstractMonoid end

# ------------------------------------------------------------------------------------
# Required functions

"""
    natoms(G::AbstractGerm) -> Int

The number `N` of atoms of `G`. The atoms are the simples `1:N`.
"""
function natoms end

"""
    nsimples(G::AbstractGerm) -> Int

The number of simples of `G`, the identity included. This number is a constant. It is defined by every tabulated germ, whose simples are then `0:nsimples(G)-1`, and by a lazy germ only if the germ knows the number without enumerating the simples.
"""
function nsimples end

# The numbers of the simples of a germ whose simples are all numbered (a tabulated germ),
# `Int32(0):Int32(nsimples(G) - 1)`: the internal method from which `simples(G)` and
# `elements(G)` are derived, and through which the construction reads the simples. Never
# defined by a lazy germ.
function _simple_numbers end

"""
    product(G::AbstractGerm, s::Integer, t::Integer) -> Union{Int32, Nothing}

The simple `st` if the product of the simples `s` and `t` is simple, and `nothing` otherwise. Where the product is known to exist, write `product(G, s, t)::Int32`: the assertion documents the invariant and checks it.
"""
function product end

"""
    ldiv(G::AbstractGerm, s::Integer, t::Integer) -> Union{Int32, Nothing}

The right complement `s \\ t` of `s` in `t`, that is, the simple `u` such that `su = t`, if `s ≼ t`, and `nothing` otherwise.
"""
function ldiv end

"""
    left_descents(G::AbstractGerm, s::Integer) -> AbstractSet{Int}

The atoms `a` such that `a ≼ s`. The result is iterated in increasing order and has a fast `in`; the germs of GarsideHomology return a [`StaticBitSet`](@ref).
"""
function left_descents end

"""
    lcm_atom(G::AbstractGerm, ω::Integer, a::Integer) -> NamedTuple{(:lcm, :complement), NTuple{2, Int32}}

The right lcm `ω ∨ a` of the simple `ω` and the atom `a`, together with the right complement `τ = ω \\ (ω ∨ a)`, that is, the simple `τ` such that `ωτ = ω ∨ a`. The result is the named tuple `(lcm = ω ∨ a, complement = τ)`, which destructures as `ωa, τ = lcm_atom(G, ω, a)`. If `a ≼ ω`, the result is `(lcm = ω, complement = 0)`.
"""
function lcm_atom end

"""
    opposite(G::AbstractGerm)

The germ of the opposite monoid `Mᵒᵖ`, which has the same atoms, in the same order, and the reversed product `s ·ᵒᵖ t = ts`. The left divisibility of `Mᵒᵖ` is the right divisibility of `M`, and its lcms are the left lcms of `M`. A tabulated opposite numbers its simples in its own canonical order, so that `opposite(opposite(G))` numbers them as `G` does.
"""
function opposite end

# ------------------------------------------------------------------------------------
# Derived functions

"""
    one(G::AbstractGerm) -> Int32

The identity of `G`, which is the simple `0`.
"""
Base.one(::AbstractGerm) = Int32(0)

"""
    isone(G::AbstractGerm, s::Integer) -> Bool

Whether the simple `s` is the identity.
"""
Base.isone(::AbstractGerm, s::Integer) = iszero(s)

"""
    atoms(G::AbstractGerm)

The atoms of `G`, that is, `Int32(1):Int32(natoms(G))`.
"""
atoms(G::AbstractGerm) = Int32(1):Int32(natoms(G))

"""
    is_left_divisor(G::AbstractGerm, s::Integer, t::Integer) -> Bool

Whether `s ≼ t`. The generic method answers without a division when `s` is the identity, when `s` is equal to `t`, when `s` is an atom (by [`left_descents`](@ref)), or when `s` is not shorter than `t`. Otherwise, it calls [`ldiv`](@ref).
"""
function is_left_divisor(G::AbstractGerm, s::Integer, t::Integer)
    iszero(s) && return true
    s == t && return true
    s <= natoms(G) && return Int(s) in left_descents(G, t)
    length(G, s) < length(G, t) || return false
    return ldiv(G, s, t) !== nothing
end

"""
    word(G::AbstractGerm, s::Integer) -> Vector{Int32}

The smallest reduced word of the simple `s` in the atoms, for the lexicographic order. It is obtained by taking, at each step, the smallest left descent of what remains. The word of the identity is empty.
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
    isless_lex(G::AbstractGerm, s::Integer, t::Integer) -> Bool

Whether the smallest word of `s` comes strictly before that of `t` in the lexicographic order, where a prefix comes first. The generic method compares the words letter by letter, without allocating; a tabulated germ compares precomputed ranks.
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

"""
    isless_deglex(G::AbstractGerm, s::Integer, t::Integer) -> Bool

Whether `s` comes strictly before `t` in the degree-lexicographic order: the shorter simple comes first, and simples of the same length are compared by their smallest words in the lexicographic order ([`isless_lex`](@ref)). This is the canonical order of the simples of a tabulated germ, whose method compares their numbers.
"""
function isless_deglex(G::AbstractGerm, s::Integer, t::Integer)
    ℓs, ℓt = length(G, s), length(G, t)
    ℓs == ℓt || return ℓs < ℓt
    return isless_lex(G, s, t)
end

# Whether concurrent queries are safe (the contract, Mutability): `true` only for a germ
# documented as immutable; a lazy germ may update its memo during a query. The tests call
# the methods of this trait, but the compiler replaces each call by its constant value, so
# that the coverage never sees their lines.
_thread_safe(::AbstractGerm) = false  # COV_EXCL_LINE

# ------------------------------------------------------------------------------------
# Entry point

"""
    germ(M::BraidMonoid)
    germ(M::DualBraidMonoid)
    germ(W::FiniteCoxeterGroup)
    germ(W::PermRootGroup)

The germ of the monoid `M`, of the type best suited to its size. The methods are added by the package extensions. With Chevie loaded, `M` may be a `BraidMonoid` of a finite Coxeter group, or a `DualBraidMonoid` of a Coxeter group or of a well-generated irreducible complex reflection group. The result depends on the number of simples, computed from the degrees:
- if this number is at most 5·10⁵, the result is a [`TableGerm`](@ref);
- above, a braid monoid gets a lazy [`CoxeterGerm`](@ref), and a dual braid monoid gives an `ArgumentError`, since there is no lazy germ for these monoids.

`TableGerm(M)` forces the tabulation, and `CoxeterGerm(cartan(W))` gives the lazy germ.

For a finite Coxeter group `W` of Chevie, `germ(W)` is the germ of its classical braid monoid, `germ(BraidMonoid(W))`. The finite Coxeter groups of Chevie are those built by `coxgroup` and `coxsym`, and the groups that `crg` builds as Coxeter groups, such as `crg(23)`, which is `coxgroup(:H, 3)`. For another reflection group of Chevie (a `PermRootGroup`, such as `crg(4)`), `germ(W)` throws an `ArgumentError`: the dual braid monoid is given by [`dual_germ`](@ref).
"""
function germ end

"""
    dual_germ(W::ReflectionGroup; c = nothing)

The germ of the dual braid monoid of the well-generated irreducible reflection group `W` for the Coxeter element `c`, given as a word in the generators of `W`. The germ depends on `c` only up to isomorphism.

By default (`c = nothing`), `c` is chosen as follows:
- `c` is the product of the generators in order, `1:n`, if this product is a regular element of order h (Springer's criterion: its eigenspace for ``e^{2i\\pi/h}`` is a line outside the reflecting hyperplanes);
- otherwise, `c` is the product of the generators in the first order, for the lexicographic order, that gives such an element, and a message names this order (for G₂₉: `[1, 2, 4, 3]`);
- if there is no such order, an `ArgumentError` is thrown.

The monoid itself is not built. For a real group (2 among the degrees), the simples are found by the criterion of Brady and Watt; for a complex group, they are found by the Hurwitz orbit of a reduced decomposition of `c` ([`hurwitz_orbit`](@ref)). For a large complex group, building `DualBraidMonoid(W)` itself with Chevie takes minutes (about 4 for G₃₄), whereas `dual_germ(W)` takes a fraction of a second.

The methods are added by the package extensions: with Chevie loaded, `W` may be a finite Coxeter group (`coxgroup`) or a complex reflection group (`crg`). The same threshold as for [`germ`](@ref) applies.
"""
function dual_germ end
