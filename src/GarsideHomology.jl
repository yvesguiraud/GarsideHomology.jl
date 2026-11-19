"""
    GarsideHomology

Small free resolutions of Garside monoids and their homology.

GarsideHomology computes the reduced lcm resolution of a Garside monoid, a variant of the resolution of Dehornoy and Lafont with fewer cells, and the homology of the monoid, which is that of its enveloping group, with coefficients in a module of rank one (for now, the trivial module and the sign module over the integers).

A computation goes through the following steps:
- a germ gives the simples of the monoid and their partial product ([`AbstractGerm`](@ref)); with Chevie loaded, [`germ`](@ref) and [`dual_germ`](@ref) build the germs of the braid monoids and of the dual braid monoids of reflection groups;
- an ordering ([`LCMOrdering`](@ref)) fixes the cells of the resolution and their matching;
- [`lcm_resolution`](@ref) builds the resolution, a chain complex of free bimodules whose bases are made of cells ([`Cell`](@ref)), and [`certify`](@ref) checks it;
- [`homology`](@ref) computes the homology with coefficients in a module ([`TrivialModule`](@ref), [`SignModule`](@ref)), or in a ring for its trivial module, as in `homology(R, ℤ)` and `homology(R, ℤε)`.

The elements of a germ are written as words in the atoms ([`GermElement`](@ref)), and the cells and the elements of the modules in bar form, with the settings of [`Notation`](@ref). The generic layer of free modules and chain complexes, `GarsideHomology.Complexes`, does not depend on the rest of the package; its exported names are also exported by GarsideHomology.
"""
module GarsideHomology

# Generic layer of free modules and chain complexes. It must not depend on anything
# else in GarsideHomology, so that it can later become a package of its own.
include("complexes/Complexes.jl")
using .Complexes
using .Complexes: rank, basis, MonoidAlgebra, AbstractMonoid, AbstractModule, AbstractChainComplex, GradedBasis,
                  BasedFreeModuleElement, Differential
using SparseArrays: nnz, nzrange, rowvals, nonzeros
import AbstractAlgebra
using AbstractAlgebra: ZZ

"""
    ℤ

The ring of integers of AbstractAlgebra, that of the coefficients of the resolutions, under a name of its own; it is the same object as `GarsideHomology.ZZ`. A ring stands for its trivial module, as in `homology(R, ℤ)` and `ℤ ⊗ R ⊗ ℤ`. The name is typed `\\bbZ` then Tab, in the Julia REPL and in the editors that know Julia (VS Code, Jupyter, Pluto, the Julia modes of Emacs and Vim).

`ZZ` is public but not exported: Nemo, Hecke and Oscar export a `ZZ` of their own, which is another object, and exporting ours would make `ZZ` ambiguous, hence unusable, for their users. To write `ZZ`, import it with `using GarsideHomology: ZZ`, rather than writing `GarsideHomology.ZZ` at each use.
"""
const ℤ = AbstractAlgebra.ZZ
import SHA
import CodecZlib
using CodecZlib: GzipCompressorStream

# Germs: numbering of the simples, partial product, divisions, descents, lcms.
include("germs/bitsets.jl")
include("germs/interface.jl")
include("germs/table.jl")
include("germs/hurwitz.jl")
include("germs/coxeter.jl")
include("germs/elements.jl")
include("germs/verify.jl")
include("germs/isomorphism.jl")

# Orderings (the preorders on the lcms that fix the cells and their matching) and
# the tree of the cells.
include("cells/orderings.jl")
include("cells/tree.jl")
include("cells/cell.jl")
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
export labels, differential, extend_scalars, torsion, is_square_zero
export homology, collapsible_pairs
export base_ring, nrows, ncols
export Notation, set_notation!, notation
# Public, not exported: Chevie exports its own `rank` and `basis` (see Complexes.jl),
# Oscar a type MonoidAlgebra of its own; GradedBasis, BasedFreeModuleElement and
# Differential, the types of basis(C), of M(x) and of differential(C), are obtained without
# being named.
public rank, basis, MonoidAlgebra, GradedBasis, BasedFreeModuleElement, Differential

# Interface of germs. Exported: the names an ordinary user types. Public only: the
# vocabulary of germs and of the construction, and the names exported elsewhere
# (`product` and `atoms` by Oscar, `word` by Chevie and Oscar).
export natoms, nsimples, simples, opposite, germ, dual_germ
public AbstractGerm, StaticBitSet
# The elements of a germ for the user: public, since one gets them by calling the
# germ (G(1, 2, 1), G("121")) or from simples(G); elements is exported by PermGroups (hence
# Chevie), FiniteFields, Hecke and Oscar (with Chevie, the extension adds a method to its
# elements).
public GermElement, elements
public product, ldiv, left_descents, lcm_atom, atoms, is_left_divisor, word, isless_lex, isless_deglex
public TableGerm, tabulate, hurwitz_orbit, CoxeterGerm, nnumbered, verify, count_multichains, GermFingerprint, fingerprint
# Tool of the developer, public only (is_isomorphic and isomorphism are exported by
# AbstractAlgebra, Hecke and Oscar).
public germ_isomorphism

# Orderings. Exported: what a user types to choose a resolution. Public only: the
# vocabulary for writing a new preorder, and the tools of the developer.
export LCMOrdering, ComplementDeglex, TrivialPreorder, FunctionPreorder
public AbstractLCMPreorder, precedes, extends_divisibility

# Tree of the cells: a tool of the construction, public but not exported, except ncells,
# which also counts the cells of a resolution. Accessors with the prefix cell_, since
# weight and degree are exported by Oscar, AbstractAlgebra, Chevie. The cells for the user:
# Cell, exported, obtained from the bases or by Cell(P, "2|12") (no Cell is exported by
# Oscar, Hecke, AbstractAlgebra or Chevie); cell_weight, cell_degree and atom_sequence,
# exported, since they also read a Cell.
export Cell, ncells, cell_weight, cell_degree, atom_sequence
public CellTree, cell_range, cell_letter, cell_atom
public cell_children, admissible_atoms, bar_form

# Projection onto the cells (provisional names, before the theory of the article).
public ChainTree, BimoduleTerm, LCMProjectionData, lcm_splitting, lcm_projection!, lcm_boundary

# The lcm resolution, its snapshots and its certificate.
export lcm_resolution, certify
public LCMResolution, CellBasis, cell_tree, write_snapshot
public ResolutionCertificate, is_correct, is_minimal, is_linear

# Modules of coefficients. Exported: the modules a user types, and ℤ, which stands for the
# trivial module. Public only: the abstract type, and the tensor products (tensor_product is
# exported by Hecke and Oscar, ⊗ by Hecke, Oscar and Chevie): write `using GarsideHomology: ⊗`
# for the notation, which wins over the exports of those packages; with Chevie loaded, its ⊗
# also tensors our complexes (package extension).
export TrivialModule, SignModule, ℤ, ℤε, ℤϵ
public AbstractMonoidModule, tensor_product, ⊗
# The ring of integers of AbstractAlgebra, the same object as ℤ: public, not exported, since
# Nemo, Hecke and Oscar export a ZZ of their own, another object.
public ZZ

# ------------------------------------------------------------------------------------
# Precompilation of the first calls: tabulated germs with at most 64 atoms
# and with more (one or two words of 64 bits for the sets of atoms, the second case being
# that of the large groups, B*(E₈) and G₃₄ to G₃₇), and a lazy germ; their resolutions for
# the two orderings, their certificates and their homology with trivial and sign
# coefficients, so that the first computation of a session does not pay for compiling the
# whole pipeline. Run at precompilation by PrecompileTools.
using PrecompileTools: @setup_workload, @compile_workload

# COV_EXCL_START: run only while the package is precompiled, never in the tests; a failure
# makes the precompilation, hence every task of the CI, fail.

# The germ of the dual braid monoid of I₂(m), tabulated from its right covers: rows, the
# atoms r₀, …, r_{m-1} (the reflections, simples 2 to m + 1); columns, the simples 1, the
# reflections and the Coxeter element c (simple m + 2); 1 · rₖ = rₖ and rₖ · rₖ₊₁ = c.
function _dual_dihedral_germ(m::Integer)
    covers = zeros(Int, m, m + 2)
    for k in 0:(m - 1)
        covers[k + 1, 1] = k + 2
        covers[mod(k + 1, m) + 1, k + 2] = m + 2
    end
    return TableGerm(covers, 1)
end

function _precompile_workload()
    # B*(A₂) and B*(I₂(65)), tabulated (65 atoms take two words), and B⁺(A₃), lazy.
    germs = (_dual_dihedral_germ(3), _dual_dihedral_germ(65), CoxeterGerm([2 -1 0; -1 2 -1; 0 -1 2]))
    for G in germs, ordering in (LCMOrdering(), LCMOrdering(TrivialPreorder()))
        R = lcm_resolution(G, ordering)
        certify(R)
        homology(R, ℤ)
        homology(R, SignModule())
        ncells(R)
    end
    # The elements, their displays and their readings, on B*(A₂) and B*(I₂(65)).
    for G in germs[1:2]
        R = lcm_resolution(G)
        d = differential(R)
        x = basis(R[last(range(R))])[1]
        v = d(x)
        R(string(v)), R(string(x)), Cell(R, string(x)), G(string(G(1))), extend_scalars(R)(string(x))
        sprint(show, MIME"text/plain"(), R; context = :limit => true)
        sprint(show, MIME"text/plain"(), R)
        sprint(show, MIME"text/plain"(), basis(R); context = :limit => true)
        sprint(show, MIME"text/plain"(), v; context = :limit => true)
        sprint(show, MIME"text/plain"(), d)
    end
    return nothing
end

@setup_workload begin
    @compile_workload _precompile_workload()
end
# COV_EXCL_STOP

end # module GarsideHomology
