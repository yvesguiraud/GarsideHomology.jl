"""
    Gauss

Small free resolutions of Garside monoids and their homology.

Gauss computes the reduced lcm resolution of a Garside monoid, a variant of the
Dehornoy–Lafont resolution with far fewer cells, and the homology of the
corresponding groups with various coefficients.

This version is a package skeleton: the layers (germs, orderings, cell tree,
projection, resolution and certificates, modules, complexes) are added step by step.
Present so far: the generic layer of complexes, `Gauss.Complexes`, whose public names
are re-exported, the interface of germs ([`AbstractGerm`](@ref)), tabulated and lazy
germs ([`TableGerm`](@ref), [`CoxeterGerm`](@ref)), the orderings
([`LCMOrdering`](@ref)), the tree of the cells ([`CellTree`](@ref)), the projection
onto the cells ([`LCMProjectionData`](@ref)) and the lcm resolution
([`lcm_resolution`](@ref)), with its snapshots ([`write_snapshot`](@ref)).
"""
module Gauss

# Generic layer of free modules and chain complexes. It must not depend on anything
# else in Gauss, so that it can later become a package of its own.
include("complexes/Complexes.jl")
using .Complexes
using .Complexes: rank, AbstractMonoid
import SHA
import CodecZlib
using CodecZlib: GzipCompressorStream

# Germs: numbering of the simples, partial product, divisions, descents, lcms.
include("germs/bitsets.jl")
include("germs/interface.jl")
include("germs/table.jl")
include("germs/hurwitz.jl")
include("germs/coxeter.jl")
include("germs/verify.jl")

# Orderings (the preorders on the lcms that fix the cells and their matching) and
# the tree of the cells.
include("cells/orderings.jl")
include("cells/tree.jl")
include("cells/chains.jl")
include("cells/projection.jl")

# The lcm resolution, a chain complex of free bimodules.
include("resolution.jl")

# Lines of the snapshots (format version 1).
include("snapshots.jl")

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
public product, ldiv, left_descents, lcm_atom, atoms, is_left_divisor, word, isless_lex, isless_deglex
public TableGerm, tabulate, hurwitz_orbit, CoxeterGerm, nnumbered, verify, count_multichains, GermFingerprint, fingerprint

# Orderings. Exported: what a user types to choose a resolution. Public only: the
# vocabulary for writing a new preorder, and the tools of the developer.
export LCMOrdering, ComplementDeglex, TrivialPreorder, FunctionPreorder
public AbstractLCMPreorder, precedes, extends_divisibility

# Tree of the cells: a tool of the construction, public but not exported. Accessors with
# the prefix cell_, since weight and degree are exported by Oscar, AbstractAlgebra, Chevie.
public CellTree, ncells, cell_range, cell_weight, cell_letter, cell_atom, cell_degree
public cell_children, admissible_atoms, bar_form, atom_sequence

# Projection onto the cells (provisional names, before the theory of the article).
public ChainTree, BimoduleTerm, LCMProjectionData, lcm_splitting, lcm_projection!, lcm_boundary

# The lcm resolution and its snapshots.
export lcm_resolution
public CellBasis, cell_tree, write_snapshot

end # module Gauss
