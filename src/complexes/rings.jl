# Rings of coefficients.
#
# Base rings are AbstractAlgebra rings (AbstractAlgebra.ZZ for now). A monoid algebra is
# described by a `MonoidAlgebra`, which keeps the monoid as an opaque object: this layer
# never multiplies elements of the monoid, it only needs to recognise the identity.

"""
    AbstractMonoid

Supertype of the objects that describe a monoid for a [`MonoidAlgebra`](@ref), such as the germs of GarsideHomology or a table of simples read from a file. Such an object fixes the numbering of the elements that serve as labels, `0` being the identity. This layer requires no method of it.
"""
abstract type AbstractMonoid end

"""
    MonoidAlgebra(base::AbstractAlgebra.Ring, monoid::AbstractMonoid, side::Symbol = :both)

The descriptor of a monoid algebra over the ring `base`, for the monoid `M` described by `monoid` ([`AbstractMonoid`](@ref)):
- with `side = :left` or `:right`, the monoid algebra `base[M]`, over which the modules are left or right modules;
- with `side = :both`, the algebra `base[M × Mᵒᵖ]`, over which the bimodules are left modules.

The elements of `M` are encoded by nonnegative integers in the numbering fixed by `monoid`, `0` being the identity; for `side = :both`, an element `(σ, τ)` of `M × Mᵒᵖ` is encoded by a pair of such integers. These codes are the labels of a [`LabelledSparseMatrix`](@ref) (see [`label_type`](@ref)). `base_ring(A)` is the ring of coefficients `base`.
"""
struct MonoidAlgebra{R <: AbstractAlgebra.Ring, M <: AbstractMonoid}
    base::R
    monoid::M
    side::Symbol
    function MonoidAlgebra(base::R, monoid::M, side::Symbol = :both) where
                {R <: AbstractAlgebra.Ring, M <: AbstractMonoid}
        side in (:left, :right, :both) ||
            throw(ArgumentError("side must be :left, :right or :both, got :$side"))
        return new{R, M}(base, monoid, side)
    end
end

# Rings over which this layer builds modules and complexes.
const _RingDescriptor = Union{AbstractAlgebra.Ring, MonoidAlgebra}

AbstractAlgebra.base_ring(A::MonoidAlgebra) = A.base

"""
    label_type(A::MonoidAlgebra)

The type of the labels of the matrices over `A`: `Int32` for a one-sided monoid algebra, and `NTuple{2, Int32}` for an algebra of bimodules.
"""
label_type(A::MonoidAlgebra) = A.side === :both ? NTuple{2, Int32} : Int32

"""
    identity_label(L::Type)

The label of the identity of the monoid, for labels of type `L`: `zero(L)` for integer labels, and `(0, 0)` for pairs.
"""
identity_label(::Type{L}) where {L <: Integer} = zero(L)
identity_label(::Type{NTuple{2, T}}) where {T <: Integer} = (zero(T), zero(T))

Base.:(==)(A::MonoidAlgebra, B::MonoidAlgebra) =
    A.side === B.side && A.base == B.base && (A.monoid === B.monoid || A.monoid == B.monoid)
Base.hash(A::MonoidAlgebra, h::UInt) = hash(A.side, hash(A.base, hash(:MonoidAlgebra, h)))

function Base.show(io::IO, A::MonoidAlgebra)
    kind = A.side === :both ? "bimodules" : "$(A.side) modules"
    print(io, "MonoidAlgebra over ", A.base, " (", kind, ", monoid of type ",
          nameof(typeof(A.monoid)), ")")
end
