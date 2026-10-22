# Cell schemes: the preorders on the lcms that fix the cells of the resolution and their
# matching, read by the construction of the cell tree (J1.d, task d2).

"""
    AbstractLcmPreorder

Supertype of the preorders of a [`CellScheme`](@ref): for each simple `ω`, a total
preorder on the right lcms `ω ∨ a` of `ω` with the atoms, read in the interval `[ω, Δ]`,
that is through their complements `ω \\ (ω ∨ a)`.

A preorder is given by [`precedes`](@ref)`(P, G, ω, x, y)`, where `x` and `y` are results
of [`lcm_atom`](@ref)`(G, ω, a)`: whether `x` comes strictly before `y`. It must be a
strict weak order on these results (irreflexive and transitive, incomparability being
transitive too). Two results neither of which precedes the other are tied: the
construction then keeps the inherited order.

The preorder is *compatible with divisibility* at `ω` when, for atoms `a` and `b`, a
strict divisibility `ω ∨ a ≺ ω ∨ b` implies that `ω ∨ a` precedes `ω ∨ b`. With a
compatible preorder, the construction keeps exactly one atom for each minimal lcm (the
reduced construction); see [`extends_divisibility`](@ref).

Methods of Gauss: [`ComplementDeglex`](@ref) (the default), [`TrivialPreorder`](@ref)
and [`FunctionPreorder`](@ref). A new preorder is a subtype with a method of `precedes`.
"""
abstract type AbstractLcmPreorder end

"""
    precedes(P::AbstractLcmPreorder, G::AbstractGerm, ω, x, y) -> Bool

Whether `x` comes strictly before `y` for the preorder `P` at the simple `ω`, where `x`
and `y` are results `(lcm = ω ∨ a, complement = ω \\ (ω ∨ a))` of
[`lcm_atom`](@ref)`(G, ω, a)` for atoms `a` that do not divide `ω`.
"""
function precedes end

"""
    ComplementDeglex() <: AbstractLcmPreorder

The degree-lexicographic order of the complements: `x` precedes `y` when
[`isless_deglex`](@ref)`(G, x.complement, y.complement)`, the shorter complement first,
then the smallest word of the complement in the lexicographic order. The default
preorder of a [`CellScheme`](@ref).

Lengths add, `ℓ(ω ∨ a) = ℓ(ω) + ℓ(ω \\ (ω ∨ a))`, so that the preorder is compatible with
divisibility; two atoms are tied only if they have the same lcm, so that the inherited
order only chooses which atom represents an lcm. It gives the reduced resolution, Gauss_a
of the preliminary implementation (`filter = true, lex = true`), on the cases of J0.
Comparing the lcms themselves in the degree-lexicographic order would not do: the
smallest word of `ω ∨ a` often starts with `a`, so that this order follows that of the
atoms rather than the shape of `[ω, Δ]`.
"""
struct ComplementDeglex <: AbstractLcmPreorder end

precedes(::ComplementDeglex, G::AbstractGerm, ω::Integer, x, y) =
    isless_deglex(G, x.complement, y.complement)

"""
    TrivialPreorder() <: AbstractLcmPreorder

The preorder in which all the lcms are tied: the construction keeps the order of the
root everywhere. With the order `p` of the atoms at the root,
`CellScheme(TrivialPreorder(); root_order = p)` is the order resolution of Dehornoy and
Lafont for `p`, Gauss_d of the preliminary implementation (`filter = false,
lex = false`). It is compatible with divisibility only where no two lcms `ω ∨ a`, `ω ∨ b`
are strictly comparable.
"""
struct TrivialPreorder <: AbstractLcmPreorder end

precedes(::TrivialPreorder, G::AbstractGerm, ω::Integer, x, y) = false

"""
    FunctionPreorder(lt) <: AbstractLcmPreorder

The preorder given by a function: `x` precedes `y` when `lt(G, ω, x, y)`, which must
return a `Bool` and be a strict weak order (see [`AbstractLcmPreorder`](@ref)). For the
exploration of other preorders; the lexicographic order of the complements, for
instance, is `FunctionPreorder((G, ω, x, y) -> Gauss.isless_lex(G, x.complement,
y.complement))`. A preorder that is not compatible with divisibility gives a resolution
that is not reduced; the result is to be checked (∂² = 0, Euler characteristic of each
weight).
"""
struct FunctionPreorder{F} <: AbstractLcmPreorder
    lt::F
end

precedes(P::FunctionPreorder, G::AbstractGerm, ω::Integer, x, y) = P.lt(G, ω, x, y)::Bool

# ------------------------------------------------------------------------------------
# Cell schemes

"""
    CellScheme(preorder = ComplementDeglex(); root_order = nothing)

The data that fix the cells of the resolution of a germ and their matching: a preorder
on the lcms ([`AbstractLcmPreorder`](@ref)) and the order of the atoms at the root, a
permutation given as the list of the atoms in increasing order (`nothing`: the order of
the germ, `1:N`).

The construction of the cells reads them as follows. At a node of weight `ω` whose list
`I` of admissible atoms is in the inherited order (the order of the root, at the root),
it sorts `I` by the preorder at `ω`, ties keeping the inherited order; it keeps an atom
`j` if no atom placed before `j` divides `ω ∨ aⱼ`; the `k`-th atom kept gives the child
of weight `ω ∨ a` whose last letter is the complement `ω \\ (ω ∨ a)` and whose list is
formed of the `k − 1` atoms kept before it, in that order. The same order serves the
matching of the cells.

- `CellScheme()`: the reduced resolution ([`ComplementDeglex`](@ref)), the default.
- `CellScheme(TrivialPreorder(); root_order = p)`: the order resolution of Dehornoy and
  Lafont for the order `p` of the atoms ([`TrivialPreorder`](@ref)).

The order of the root acts only through the ties of the preorder. With
`ComplementDeglex`, ties join atoms with the same lcm, and the order of the root only
changes which atom represents it.
"""
struct CellScheme{P <: AbstractLcmPreorder}
    preorder::P
    root_order::Union{Nothing, Vector{Int32}}
    function CellScheme{P}(preorder::P, root_order) where {P <: AbstractLcmPreorder}
        return new{P}(preorder, root_order === nothing ? nothing : _permutation(root_order))
    end
end

CellScheme(preorder::P = ComplementDeglex(); root_order = nothing) where {P <: AbstractLcmPreorder} =
    CellScheme{P}(preorder, root_order)

# A permutation of 1:n given as a list, checked and copied as Int32.
function _permutation(p)
    q = Int32[]
    for x in p
        x isa Integer || throw(ArgumentError("CellScheme: root_order contains $x, not an integer"))
        push!(q, Int32(x))
    end
    sort(q) == 1:length(q) ||
        throw(ArgumentError("CellScheme: root_order = $(Int.(q)) is not a permutation of 1:$(length(q))"))
    return q
end

# The atoms of G in the order of the root of the scheme S.
function _root_order(S::CellScheme, G::AbstractGerm)
    N = natoms(G)
    S.root_order === nothing && return collect(Int32(1):Int32(N))
    length(S.root_order) == N ||
        throw(ArgumentError("CellScheme: root_order has $(length(S.root_order)) atoms, the germ $N"))
    return copy(S.root_order)
end

Base.:(==)(S::CellScheme, T::CellScheme) = S.preorder == T.preorder && S.root_order == T.root_order
Base.hash(S::CellScheme, h::UInt) = hash(S.root_order, hash(S.preorder, hash(:CellScheme, h)))

function Base.show(io::IO, S::CellScheme)
    print(io, "CellScheme(")
    show(io, S.preorder)
    S.root_order === nothing || print(io, "; root_order = ", Int.(S.root_order))
    print(io, ")")
end

# ------------------------------------------------------------------------------------
# Compatibility with divisibility (a tool of the developer)

"""
    extends_divisibility(P::AbstractLcmPreorder, G::AbstractGerm, ωs = simples(G)) -> Bool

Whether the preorder `P` is compatible with divisibility at each simple of `ωs`: for all
atoms `a` and `b` that do not divide `ω`, if `ω ∨ a` strictly divides `ω ∨ b`, then
`ω ∨ a` precedes `ω ∨ b`. A scheme whose preorder is compatible keeps exactly one atom
for each minimal lcm. The test does not depend on the tree of the construction: an
incompatible preorder may still give a reduced resolution, when the pairs that break the
condition never meet in the same list.

By brute force, in `N²` comparisons for each simple (`N` atoms): a tool of the developer.
For a lazy germ, which does not list its simples, give them as `ωs`.
"""
function extends_divisibility(P::AbstractLcmPreorder, G::AbstractGerm, ωs = simples(G))
    N = natoms(G)
    for ω in ωs
        D = left_descents(G, ω)
        for a in 1:N
            a in D && continue
            x = lcm_atom(G, ω, a)
            for b in 1:N
                (b == a || b in D) && continue
                y = lcm_atom(G, ω, b)
                # ω ∨ a ≼ ω ∨ b if and only if a ≼ ω ∨ b, since ω ≼ ω ∨ b.
                x.lcm != y.lcm && a in left_descents(G, y.lcm) || continue
                precedes(P, G, ω, x, y) || return false
            end
        end
    end
    return true
end
