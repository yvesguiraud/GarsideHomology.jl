# Rings of coefficients.
#
# Base rings are AbstractAlgebra rings (AbstractAlgebra.ZZ for now). A monoid algebra is
# described by a `MonoidAlgebra`, which keeps the monoid as an opaque object: this layer
# never multiplies elements of the monoid, it only needs to recognise the identity.

"""
    AbstractMonoid

Supertype of the objects that describe a monoid for a [`MonoidAlgebra`](@ref): the germs
of the Gauss layer, a table of simples read from a file, etc. Such an object fixes the
numbering of the elements used as labels (1 for the identity). No method is required
by this layer.
"""
abstract type AbstractMonoid end

"""
    MonoidAlgebra(base, monoid, side = :both)

Descriptor of the monoid algebra `base[M]` of a monoid `M` (`side = :left` or `:right`),
or of the algebra `base[M × Mᵒᵖ]` over which bimodules are left modules (`side = :both`).
`base` is an AbstractAlgebra ring and `monoid` an [`AbstractMonoid`](@ref).

Elements of `M` are encoded by positive integers in the numbering fixed by `monoid`,
`1` being the identity; for `side = :both`, an element `(σ, τ)` of `M × Mᵒᵖ` is encoded
by a pair of such integers. These codes are the labels of a
[`LabelledSparseMatrix`](@ref) (see [`label_type`](@ref)).

`base_ring(A)` is the ring of coefficients `base`.
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

Type of the labels of the matrices over `A`: `Int32` for a one-sided monoid algebra,
`NTuple{2, Int32}` for an algebra of bimodules.
"""
label_type(A::MonoidAlgebra) = A.side === :both ? NTuple{2, Int32} : Int32

"""
    identity_label(L)

Label of the identity of the monoid: `one(L)` for integer labels, `(1, 1)` for pairs.
"""
identity_label(::Type{L}) where {L <: Integer} = one(L)
identity_label(::Type{NTuple{2, T}}) where {T <: Integer} = (one(T), one(T))

Base.:(==)(A::MonoidAlgebra, B::MonoidAlgebra) =
    A.side === B.side && A.base == B.base && (A.monoid === B.monoid || A.monoid == B.monoid)
Base.hash(A::MonoidAlgebra, h::UInt) = hash(A.side, hash(A.base, hash(:MonoidAlgebra, h)))

function Base.show(io::IO, A::MonoidAlgebra)
    kind = A.side === :both ? "bimodules" : "$(A.side) modules"
    print(io, "MonoidAlgebra over ", A.base, " (", kind, ", monoid of type ",
          nameof(typeof(A.monoid)), ")")
end
