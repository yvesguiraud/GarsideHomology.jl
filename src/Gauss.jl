"""
    Gauss

Small free resolutions of Garside monoids and their homology.

Gauss computes the reduced lcm resolution of a Garside monoid, a variant of the
Dehornoy–Lafont resolution with far fewer cells, and the homology of the
corresponding groups with various coefficients.

This version is a package skeleton: the layers (germs, cell schemes, cell tree,
projection, resolution and certificates, modules, complexes) are added step by step.
Present so far: the generic layer of complexes, `Gauss.Complexes`, whose public names
are re-exported, the interface of germs ([`AbstractGerm`](@ref)) and tabulated germs
([`TableGerm`](@ref)).
"""
module Gauss

# Generic layer of free modules and chain complexes. It must not depend on anything
# else in Gauss, so that it can later become a package of its own.
include("complexes/Complexes.jl")
using .Complexes
using .Complexes: rank, AbstractMonoid

# Germs: numbering of the simples, partial product, divisions, descents, lcms.
include("germs/bitsets.jl")
include("germs/interface.jl")
include("germs/table.jl")
include("germs/hurwitz.jl")

# Public interface of the complexes layer (kept equal to its export list by the tests).
export MonoidAlgebra, BasedFreeModule, LabelledSparseMatrix, ChainComplex
export FinitelyGeneratedModule, FinitelyGeneratedAbelianGroup
export FinitelyGeneratedGradedModule, FinitelyGeneratedGradedAbelianGroup
export labelled_sparse, cells, labels, differential, specialise, torsion, is_square_zero
export homology, collapsible_pairs
export base_ring, nrows, ncols
# Public, not exported: Chevie exports its own `rank` (see Complexes.jl).
public rank

# Interface of germs. Exported: the names an ordinary user types. Public only: the
# vocabulary of germs and of the construction, and the names exported elsewhere
# (`product` and `atoms` by Oscar, `word` by Chevie and Oscar).
export natoms, nsimples, simples, opposite, germ, dual_germ
public AbstractGerm, StaticBitSet
public product, ldiv, left_descents, lcm_atom, atoms, is_left_divisor, word, isless_lex
public TableGerm, tabulate, hurwitz_orbit

end # module Gauss
