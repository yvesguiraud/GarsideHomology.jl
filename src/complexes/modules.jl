# Free modules with a chosen basis.

"""
    AbstractModule{R}

Supertype of the modules of this layer, over a ring described by an object of type `R`. The interface is made of `base_ring(M)`, the ring, and `rank(M)`, the number of elements of a basis.
"""
abstract type AbstractModule{R} end

"""
    BasedFreeModule(ring::Union{AbstractAlgebra.Ring, MonoidAlgebra}, basis::AbstractVector)
    BasedFreeModule(ring::Union{AbstractAlgebra.Ring, MonoidAlgebra}, n::Integer)

The free module over `ring`, an AbstractAlgebra ring or a [`MonoidAlgebra`](@ref), with the chosen basis `basis`, indexed by `1:length(basis)`. The elements of `basis` describe the elements of the basis (names, cells, …). The basis is shared, never copied, between the complexes that use it, such as a complex and its extensions of scalars. `BasedFreeModule(ring, n)` takes `Base.OneTo(n)` as basis.

The accessors are `rank(M)`, the number of elements of the basis, [`basis`](@ref)`(M)`, the basis, and `base_ring(M)`, the ring.
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

The basis of `M`: the descriptors of the elements of its basis, in their order. This is a method of the function `basis` of AbstractAlgebra. The name is public but not exported, because Chevie exports a `basis` of its own: write `GarsideHomology.basis(M)`, or `using GarsideHomology: basis` when Chevie is not loaded. With Chevie loaded, `basis(M)` works without qualification, since the extension adds methods to the `basis` of Chevie.
"""
AbstractAlgebra.basis(M::BasedFreeModule) = M.basis

Base.:(==)(M::BasedFreeModule, N::BasedFreeModule) = M.ring == N.ring && M.basis == N.basis
Base.hash(M::BasedFreeModule, h::UInt) = hash(M.basis, hash(M.ring, hash(:BasedFreeModule, h)))

Base.show(io::IO, M::BasedFreeModule) =
    print(io, "BasedFreeModule of rank ", rank(M), " over ", M.ring)
