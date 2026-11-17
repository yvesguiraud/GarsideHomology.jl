# Modules of coefficients: modules over the monoid algebra of a germ, free of finite rank
# over a ring, given by the action of the atoms.

"""
    AbstractMonoidModule{R}

Supertype of the modules of coefficients. Such a module is a module over the monoid algebra `R[M]` of a monoid `M`, free of finite rank over the ring `R` (an AbstractAlgebra ring), and given by the action of the atoms; a simple acts as a word of atoms. A module does not depend on a germ: the tensor products ([`tensor_product`](@ref)) compute the action on the simples of the germ of the complex, and check, when needed, that it is well defined.

The interface is made of `base_ring(M)`, the ring `R`, and `rank(M)`, the rank over `R`. A module of rank one is given by a character χ of the monoid, and serves on both sides: `s·m = m·s = χ(s)·m`.
"""
abstract type AbstractMonoidModule{R <: AbstractAlgebra.Ring} end

"""
    TrivialModule(R::AbstractAlgebra.Ring = GarsideHomology.ZZ)

The trivial module: the ring `R`, on which every element of the monoid acts by 1. By default, `R` is the ring of integers of AbstractAlgebra (`GarsideHomology.ZZ`), that of the coefficients of the resolutions. The homology with coefficients in `TrivialModule()` is the homology of the monoid with trivial integer coefficients, which is also that of its enveloping group, since the monoid is a Garside monoid.
"""
struct TrivialModule{R} <: AbstractMonoidModule{R}
    ring::R
    TrivialModule(R::AbstractAlgebra.Ring) = new{typeof(R)}(R)
end

TrivialModule() = TrivialModule(AbstractAlgebra.ZZ)

"""
    SignModule(R::AbstractAlgebra.Ring = GarsideHomology.ZZ)

The sign module: the ring `R`, on which every atom acts by −1, so that a simple `s` acts by `(−1)^ℓ(s)`. For a braid group and the integers, this is the module ℤε in which the braid reflections act by −1 (Garnier, table 4). The module is well defined for every homogeneous germ. By default, `R` is the ring of integers of AbstractAlgebra (`GarsideHomology.ZZ`).
"""
struct SignModule{R} <: AbstractMonoidModule{R}
    ring::R
    SignModule(R::AbstractAlgebra.Ring) = new{typeof(R)}(R)
end

SignModule() = SignModule(AbstractAlgebra.ZZ)

"""
    ℤε

The sign module over the integers, `SignModule()`, on which every atom acts by −1. With it, `homology(R, ℤε)` is the homology with coefficients in ℤε (Garnier, table 4), beside `homology(R, ℤ)` for the trivial coefficients, and `ℤ ⊗ R ⊗ ℤε` is the complex that computes it.

The name is typed `\\bbZ` then Tab, then `\\varepsilon` then Tab (ε, U+03B5). The other epsilon, `\\epsilon` (ϵ, U+03F5, another character), gives [`ℤϵ`](@ref), the same object.
"""
const ℤε = SignModule()

"""
    ℤϵ

The same object as [`ℤε`](@ref), the sign module over the integers, under the name typed with `\\epsilon` (ϵ, U+03F5) rather than `\\varepsilon` (ε, U+03B5): these two characters look alike.
"""
const ℤϵ = ℤε

const _RankOneModule = Union{TrivialModule, SignModule}

AbstractAlgebra.base_ring(M::_RankOneModule) = M.ring
# The tests call this method, but the compiler replaces each call by the constant 1, so
# that the coverage never sees the line (as for precedes and _thread_safe).
AbstractAlgebra.rank(::_RankOneModule) = 1  # COV_EXCL_LINE

Base.show(io::IO, M::TrivialModule) = print(io, "TrivialModule over ", M.ring)
Base.show(io::IO, M::SignModule) = print(io, "SignModule over ", M.ring)

# The character of a module of rank one on the simples of the germ G, with integer values.
_character(::TrivialModule, G::AbstractGerm) = s -> 1
_character(::SignModule, G::AbstractGerm) = s -> isodd(length(G, s)) ? -1 : 1
