# Free modules with a chosen basis.

"""
    AbstractModule{R}

Supertype of the modules of this layer, over a ring described by an object of type `R`.
Interface: `base_ring(M)`, the ring, and `rank(M)`, the number of elements of a basis.
"""
abstract type AbstractModule{R} end

"""
    BasedFreeModule(ring, basis)
    BasedFreeModule(ring, n::Integer)

Free module over `ring` (an AbstractAlgebra ring or a [`MonoidAlgebra`](@ref)) with the
chosen basis `basis`, indexed by `1:length(basis)`. The elements of `basis` describe the
basis elements (names, atom sequences, ...); the basis is shared, never copied, between
the complexes that use it, e.g. a complex and its extensions of scalars.
`BasedFreeModule(ring, n)` takes `Base.OneTo(n)` as basis.

`rank(M)` is the number of basis elements, [`basis`](@ref)`(M)` the basis, `base_ring(M)`
the ring.
"""
struct BasedFreeModule{R <: _RingDescriptor, B <: AbstractVector} <: AbstractModule{R}
    ring::R
    basis::B
end

function BasedFreeModule(ring::_RingDescriptor, n::Integer)
    n >= 0 || throw(ArgumentError("negative rank $n"))
    return BasedFreeModule(ring, Base.OneTo(Int(n)))
end

AbstractAlgebra.base_ring(M::BasedFreeModule) = M.ring
AbstractAlgebra.rank(M::BasedFreeModule) = length(M.basis)

"""
    basis(M::BasedFreeModule)

The basis of `M`: the descriptors of its basis elements, in their order. A method of the
function `basis` of AbstractAlgebra, public but not exported, because Chevie exports a
`basis` of its own: write `GarsideHomology.basis(M)`, or `using GarsideHomology: basis`
when Chevie is not loaded.
"""
AbstractAlgebra.basis(M::BasedFreeModule) = M.basis

Base.:(==)(M::BasedFreeModule, N::BasedFreeModule) = M.ring == N.ring && M.basis == N.basis
Base.hash(M::BasedFreeModule, h::UInt) = hash(M.basis, hash(M.ring, hash(:BasedFreeModule, h)))

Base.show(io::IO, M::BasedFreeModule) =
    print(io, "BasedFreeModule of rank ", rank(M), " over ", M.ring)
