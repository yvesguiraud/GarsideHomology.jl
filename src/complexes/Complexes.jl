"""
    GarsideHomology.Complexes

Free modules with a basis, labelled sparse matrices, chain complexes and their homology.

The rings of coefficients are AbstractAlgebra rings (`AbstractAlgebra.ZZ` for now) and
monoid algebras described by [`MonoidAlgebra`](@ref). A complex over a monoid algebra
has [`LabelledSparseMatrix`](@ref) maps; a complex over ℤ has `SparseMatrixCSC` maps.
The homology of a complex over ℤ, computed by `SparseElementaryDivisors`, is a
[`FinitelyGeneratedGradedAbelianGroup`](@ref) of [`FinitelyGeneratedAbelianGroup`](@ref)s.

This submodule depends on nothing else in GarsideHomology, and its tests build no monoid,
so that it can later become a separate package.
"""
module Complexes

using SparseArrays: SparseArrays, SparseMatrixCSC, spzeros, nnz, nzrange, rowvals, nonzeros
import AbstractAlgebra
import SparseElementaryDivisors
using SparseElementaryDivisors: sparse_elementary_divisors
# The same functions as AbstractAlgebra's: methods are added for our types, so that the
# names do not clash with AbstractAlgebra, Nemo, Hecke or Oscar. `rank` and `basis` are
# public but not exported, because Chevie exports a `rank` (that of Chevie.Symbols, the
# rank of a partition) and a `basis` (that of Chevie.Algebras) of its own: with Chevie
# loaded, the package extension adds methods for our types to those functions, so that
# `rank(M)` and `basis(M)` work without qualification (J1.h); without Chevie, write
# `using GarsideHomology: rank, basis`.
using AbstractAlgebra: base_ring, rank, basis, nrows, ncols

export BasedFreeModule, LabelledSparseMatrix, ChainComplex
export FinitelyGeneratedModule, FinitelyGeneratedAbelianGroup
export FinitelyGeneratedGradedModule, FinitelyGeneratedGradedAbelianGroup
export labels, differential, extend_scalars, torsion, is_square_zero
export homology, collapsible_pairs
export base_ring, nrows, ncols
export Notation, set_notation!, notation

# MonoidAlgebra is public but not exported: Oscar exports a type of its own of that name.
public rank, basis, MonoidAlgebra, AbstractMonoid, AbstractModule, AbstractChainComplex, label_type, identity_label, check_invariants

include("arithmetic.jl")
include("notation.jl")
include("rings.jl")
include("modules.jl")
include("labelled.jl")
include("complex.jl")
include("fg_modules.jl")
include("homology.jl")
include("pairs.jl")

end # module Complexes
