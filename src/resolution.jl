# The lcm resolution of a germ, as a chain complex of free bimodules (J1.d, task d4).

"""
    CellBasis

The basis of the module of degree `n` of an lcm resolution: the cells of degree `n` of
its tree, as their numbers `cell_range(T, n)`, in the order of the construction. A
vector of `Int` that keeps the tree, so that a complex built on such bases, and its
specialisations, which share them, know their cells ([`cell_tree`](@ref),
[`atom_sequence`](@ref), [`bar_form`](@ref)).
"""
struct CellBasis{T <: CellTree} <: AbstractVector{Int}
    tree::T
    degree::Int
end

Base.size(B::CellBasis) = (ncells(B.tree, B.degree),)
Base.IndexStyle(::Type{<:CellBasis}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(B::CellBasis, i::Int)
    @boundscheck checkbounds(B, i)
    return first(cell_range(B.tree, B.degree)) + i - 1
end

Base.show(io::IO, B::CellBasis) = print(io, "CellBasis: ", length(B), " cells of degree ", B.degree)
Base.show(io::IO, ::MIME"text/plain", B::CellBasis) = show(io, B)

"""
    lcm_resolution(G::AbstractGerm, ordering::LCMOrdering = LCMOrdering())
    lcm_resolution(P::LCMProjectionData)

The lcm resolution of the germ `G` for the ordering `ordering`: the free resolution of
the monoid algebra `ℤ[M]` by `ℤ[M]`-bimodules obtained from the Garside complex by the
matching of the cells (see [`LCMOrdering`](@ref), [`CellTree`](@ref),
[`LCMProjectionData`](@ref)), as a `ChainComplex` over `MonoidAlgebra(ZZ, G, :both)`.

Its module of degree `n` has for basis the cells of degree `n` ([`CellBasis`](@ref)), and
`differential(C, n)`, from degree `n + 1` to degree `n`, is a `LabelledSparseMatrix` whose
column `j` is the boundary of the `j`-th cell of degree `n + 1`, with labels `(σ, τ)` for
the terms `coef · σ y τ`. Its specialisations ([`specialise`](@ref)) compute the homology
of the monoid, that of its group of fractions.

- `lcm_resolution(G)`: the reduced resolution ([`ComplementDeglex`](@ref)).
- `lcm_resolution(G, LCMOrdering(TrivialPreorder(); root_order = p))`: the order
  resolution of Dehornoy and Lafont for the order `p` of the atoms.

The first form drops the memo of the projection at the end. The second one reads the
boundaries from data `P` built by the caller, who keeps its memo.
"""
lcm_resolution(G::AbstractGerm, ordering::LCMOrdering = LCMOrdering()) =
    lcm_resolution(LCMProjectionData(CellTree(G, ordering)))

function lcm_resolution(P::LCMProjectionData)
    T = P.tree
    A = MonoidAlgebra(Complexes.AbstractAlgebra.ZZ, T.germ, :both)
    top = _top_degree(T)
    modules = [BasedFreeModule(A, CellBasis(T, n)) for n in 0:top]
    maps = LabelledSparseMatrix{NTuple{2, Int32}, Int}[_differential_matrix(P, n) for n in 0:(top - 1)]
    return ChainComplex(modules, maps; lo = 0)
end

# The matrix of d_n, from the cells of degree n + 1 (columns) to those of degree n (rows),
# filled column by column from the boundaries, whose normal form is the canonical storage.
function _differential_matrix(P::LCMProjectionData, n::Integer)
    T = P.tree
    rows, columns = cell_range(T, n), cell_range(T, n + 1)
    first_row = first(rows)
    colptr = Vector{Int}(undef, length(columns) + 1)
    colptr[1] = 1
    rowval = Int32[]
    label = NTuple{2, Int32}[]
    nzval = Int[]
    buffer = BimoduleTerm[]
    for (j, x) in enumerate(columns)
        _boundary!(buffer, P, x)
        for t in buffer
            push!(rowval, Int32(t.cell - first_row + 1))
            push!(label, (t.left, t.right))
            push!(nzval, Int(t.coef))
        end
        colptr[j + 1] = length(rowval) + 1
    end
    return LabelledSparseMatrix(length(rows), length(columns), colptr, rowval, label, nzval)
end

"""
    cell_tree(C::ChainComplex) -> CellTree

The tree of the cells of an lcm resolution, or of one of its specialisations, read from
the bases of its modules.
"""
function cell_tree(C::ChainComplex)
    b = cells(C[first(range(C))])
    b isa CellBasis || throw(ArgumentError("cell_tree: the complex is not an lcm resolution (its bases are not CellBasis)"))
    return b.tree
end
