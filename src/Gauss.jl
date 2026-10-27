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
([`lcm_resolution`](@ref), [`LCMResolution`](@ref)), with its snapshots
([`write_snapshot`](@ref)) and its certificate ([`certify`](@ref)), the modules of
coefficients ([`TrivialModule`](@ref), [`SignModule`](@ref)) and the tensor products by
them ([`tensor_product`](@ref), `homology(R, M)`).
"""
module Gauss

# Generic layer of free modules and chain complexes. It must not depend on anything
# else in Gauss, so that it can later become a package of its own.
include("complexes/Complexes.jl")
using .Complexes
using .Complexes: rank, MonoidAlgebra, AbstractMonoid, AbstractModule, AbstractChainComplex
using SparseArrays: nnz, nzrange, rowvals, nonzeros
import AbstractAlgebra
using AbstractAlgebra: ZZ
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
include("germs/isomorphism.jl")

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

# Certificates of the resolutions: correctness, minimality, linearity.
include("certificates.jl")

# Modules of coefficients and tensor products over the monoid algebra.
include("modules/monoid_modules.jl")
include("modules/tensor.jl")

# Public interface of the complexes layer (kept equal to its export list by the tests).
export BasedFreeModule, LabelledSparseMatrix, ChainComplex
export FinitelyGeneratedModule, FinitelyGeneratedAbelianGroup
export FinitelyGeneratedGradedModule, FinitelyGeneratedGradedAbelianGroup
export labelled_sparse, cells, labels, differential, specialise, torsion, is_square_zero
export homology, collapsible_pairs
export base_ring, nrows, ncols
# Public, not exported: Chevie exports its own `rank` (see Complexes.jl), Oscar a type
# MonoidAlgebra of its own (J1.e, M5).
public rank, MonoidAlgebra

# Interface of germs. Exported: the names an ordinary user types. Public only: the
# vocabulary of germs and of the construction, and the names exported elsewhere
# (`product` and `atoms` by Oscar, `word` by Chevie and Oscar).
export natoms, nsimples, simples, opposite, germ, dual_germ
public AbstractGerm, StaticBitSet
public product, ldiv, left_descents, lcm_atom, atoms, is_left_divisor, word, isless_lex, isless_deglex
public TableGerm, tabulate, hurwitz_orbit, CoxeterGerm, nnumbered, verify, count_multichains, GermFingerprint, fingerprint
# Tool of the developer, public only (is_isomorphic and isomorphism are exported by
# AbstractAlgebra, Hecke and Oscar).
public germ_isomorphism

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

# The lcm resolution, its snapshots and its certificate.
export lcm_resolution, certify
public LCMResolution, CellBasis, cell_tree, write_snapshot
public ResolutionCertificate, is_correct, is_minimal, is_linear

# Modules of coefficients. Exported: the modules a user types. Public only: the abstract
# type, and the tensor products (tensor_product is exported by Hecke and Oscar, ⊗ by Hecke):
# write `using Gauss: ⊗` for the notation, which wins over the exports of those packages.
export TrivialModule, SignModule
public AbstractMonoidModule, tensor_product, ⊗
# The ring of integers of AbstractAlgebra, that of the coefficients of the resolutions:
# public, not exported, since Nemo, Hecke and Oscar export a ZZ of their own, another object.
public ZZ

end # module Gauss
