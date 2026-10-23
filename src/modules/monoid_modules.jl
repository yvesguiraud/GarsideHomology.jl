# Modules of coefficients: modules over the monoid algebra of a germ, free of finite rank
# over a ring, given by the action of the atoms (J1.e, task e2).

"""
    AbstractMonoidModule{R}

Supertype of the modules of coefficients: a module over the monoid algebra `R[M]` of a
monoid `M`, free of finite rank over the ring `R` (an AbstractAlgebra ring), given by the
action of the atoms; a simple acts as a word of atoms. A module does not depend on a germ:
the tensor products ([`tensor_product`](@ref)) compute the action on the simples of the
germ of the complex (and check, when needed, that it is well defined).

Interface: `base_ring(M)`, the ring `R`, and `rank(M)`, the rank over `R`. A module of rank
one is given by a character χ of the monoid, and serves on both sides, `s·m = m·s = χ(s)·m`.
"""
abstract type AbstractMonoidModule{R <: AbstractAlgebra.Ring} end

"""
    TrivialModule(R = Gauss.ZZ)

The ring `R`, every element of the monoid acting by 1: the trivial module. The homology
with coefficients in `TrivialModule()` is the homology of the group of fractions with
trivial integer coefficients. By default, `R` is the ring of integers of AbstractAlgebra
(`Gauss.ZZ`), that of the coefficients of the resolutions.
"""
struct TrivialModule{R} <: AbstractMonoidModule{R}
    ring::R
    TrivialModule(R::AbstractAlgebra.Ring) = new{typeof(R)}(R)
end

TrivialModule() = TrivialModule(AbstractAlgebra.ZZ)

"""
    SignModule(R = Gauss.ZZ)

The ring `R`, every atom acting by −1, so that a simple `s` acts by `(−1)^ℓ(s)`: for a braid
group and the integers, the module ℤε in which the braid reflections act by −1 (Garnier,
table 4). Well defined on every homogeneous germ. By default, `R` is the ring of integers
of AbstractAlgebra (`Gauss.ZZ`).
"""
struct SignModule{R} <: AbstractMonoidModule{R}
    ring::R
    SignModule(R::AbstractAlgebra.Ring) = new{typeof(R)}(R)
end

SignModule() = SignModule(AbstractAlgebra.ZZ)

const _RankOneModule = Union{TrivialModule, SignModule}

AbstractAlgebra.base_ring(M::_RankOneModule) = M.ring
AbstractAlgebra.rank(::_RankOneModule) = 1

Base.show(io::IO, M::TrivialModule) = print(io, "TrivialModule over ", M.ring)
Base.show(io::IO, M::SignModule) = print(io, "SignModule over ", M.ring)

# The character of a module of rank one on the simples of the germ G, with integer values.
_character(::TrivialModule, G::AbstractGerm) = s -> 1
_character(::SignModule, G::AbstractGerm) = s -> isodd(length(G, s)) ? -1 : 1
