# Free modules with a basis of cells.

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
cells (names, atom sequences, ...); the basis is shared, never copied, between the
complexes that use it, e.g. a complex and its specialisations.
`BasedFreeModule(ring, n)` takes `Base.OneTo(n)` as basis.

`rank(M)` is the number of cells, [`cells`](@ref)`(M)` the basis, `base_ring(M)` the ring.
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
    cells(M::BasedFreeModule)

The basis of `M`: the descriptors of its cells, in the order of the basis.
"""
cells(M::BasedFreeModule) = M.basis

Base.:(==)(M::BasedFreeModule, N::BasedFreeModule) = M.ring == N.ring && M.basis == N.basis
Base.hash(M::BasedFreeModule, h::UInt) = hash(M.basis, hash(M.ring, hash(:BasedFreeModule, h)))

Base.show(io::IO, M::BasedFreeModule) =
    print(io, "BasedFreeModule of rank ", rank(M), " over ", M.ring)
