"""
    Gauss.Complexes

Free modules with a basis, labelled sparse matrices, chain complexes and their homology.

The rings of coefficients are AbstractAlgebra rings (`AbstractAlgebra.ZZ` for now) and
monoid algebras described by [`MonoidAlgebra`](@ref). A complex over a monoid algebra
has [`LabelledSparseMatrix`](@ref) maps; a complex over ℤ has `SparseMatrixCSC` maps.
The homology of a complex over ℤ, computed by `SparseElementaryDivisors`, is a
[`FinitelyGeneratedGradedAbelianGroup`](@ref) of [`FinitelyGeneratedAbelianGroup`](@ref)s.

This submodule depends on nothing else in Gauss, and its tests build no monoid,
so that it can later become a separate package.
"""
module Complexes

using SparseArrays: SparseArrays, SparseMatrixCSC, spzeros, nnz, nzrange, rowvals, nonzeros
import AbstractAlgebra
import SparseElementaryDivisors
using SparseElementaryDivisors: sparse_elementary_divisors
# The same functions as AbstractAlgebra's: methods are added for our types, so that the
# names do not clash with AbstractAlgebra, Nemo, Hecke or Oscar. `rank` is public but not
# exported, because Chevie exports a function `rank` of its own (rank of a reflection
# group): write `Gauss.rank`, or `using Gauss: rank` when Chevie is not loaded.
using AbstractAlgebra: base_ring, rank, nrows, ncols

export MonoidAlgebra, BasedFreeModule, LabelledSparseMatrix, ChainComplex
export FinitelyGeneratedModule, FinitelyGeneratedAbelianGroup
export FinitelyGeneratedGradedModule, FinitelyGeneratedGradedAbelianGroup
export labelled_sparse, cells, labels, differential, specialise, torsion, is_square_zero
export homology, collapsible_pairs
export base_ring, nrows, ncols

public rank, AbstractMonoid, AbstractModule, label_type, identity_label, check_invariants

include("arithmetic.jl")
include("rings.jl")
include("modules.jl")
include("labelled.jl")
include("complex.jl")
include("fg_modules.jl")
include("homology.jl")
include("pairs.jl")

end # module Complexes
