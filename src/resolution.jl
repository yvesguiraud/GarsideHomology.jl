# The lcm resolution of a germ, as a chain complex of free bimodules.

"""
    CellBasis

The basis of the module of degree `n` of an lcm resolution: the cells of degree `n` of its tree, as [`Cell`](@ref)s, in the order of the construction, that is, the numbers `cell_range(T, n)`. The cells are created on demand.

A `CellBasis` is a vector that keeps the tree, so that a complex built on such bases, and its extensions of scalars, which share them, know their cells ([`cell_tree`](@ref)). At the REPL, it is shown as a vector of cells in bar form, cut in the middle when it is long.
"""
struct CellBasis{T <: CellTree} <: AbstractVector{Cell{T}}
    tree::T
    degree::Int
end

Base.size(B::CellBasis) = (ncells(B.tree, B.degree),)
Base.IndexStyle(::Type{<:CellBasis}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(B::CellBasis, i::Int)
    @boundscheck checkbounds(B, i)
    return Cell(B.tree, Int32(first(cell_range(B.tree, B.degree)) + i - 1))
end

Base.showarg(io::IO, B::CellBasis, toplevel::Bool) = print(io, "CellBasis of degree ", B.degree)
Base.show(io::IO, B::CellBasis) = print(io, "CellBasis of degree ", B.degree, " with ", length(B), " cells")

# Two bases of cells of the same tree and degree are equal without comparing the cells,
# which are then the same; otherwise only two empty bases are equal, as vectors.
Base.:(==)(A::CellBasis, B::CellBasis) =
    A.tree === B.tree && A.degree == B.degree || isempty(A) && isempty(B)

# The hooks of the bases of the complexes layer for the cells: a cell is written in bar form
# between brackets, [] for the cell of degree 0; it is found from itself or from its bar
# form, and it knows its degree.

struct _CellGenerators{W <: _CellWriter}
    cells::W
    empty::String
end

Complexes._generator_writer(B::CellBasis, N::Notation) = _CellGenerators(_CellWriter(B.tree, N), N.empty_cell_in_module)

function Complexes._write_generator(io::IO, w::_CellGenerators, B::CellBasis, i::Integer)
    B.degree == 0 && return print(io, w.empty)
    print(io, "[")
    _write_cell(io, w.cells, B.tree, first(cell_range(B.tree, B.degree)) + i - 1)
    print(io, "]")
    return nothing
end

function Complexes._basis_index(B::CellBasis, x::Cell)
    x.tree === B.tree && cell_degree(x) == B.degree ||
        throw(ArgumentError("the cell $x is not one of the basis (the cells of degree $(B.degree) of its tree)"))
    return Int(x.index) - first(cell_range(B.tree, B.degree)) + 1
end

# The cell written `text`, the cell of degree 0 for an empty text (the generator of []).
_written_cell(T::CellTree, N::Notation, text::AbstractString) =
    isempty(strip(text)) ? Cell(T, Int32(1)) : _read_cell(T, N, text)

Complexes._find_in_basis(B::CellBasis, N::Notation, text::AbstractString, ::AbstractString) =
    Complexes._basis_index(B, _written_cell(B.tree, N, text))

Complexes._basis_degree(B::CellBasis) = B.degree

# A complex whose bases are cells reads its generators as bar forms (P("2|12") is [2|12]).
function Complexes._find_generator(B::CellBasis, ::AbstractChainComplex, N::Notation, text::AbstractString)
    c = _written_cell(B.tree, N, text)
    n = cell_degree(c)
    return (n, Int(c.index) - first(cell_range(B.tree, n)) + 1)
end

Complexes._generator_degree(::CellBasis, x::Cell) = cell_degree(x)

# The detailed displays (display.jl of the complexes layer): the section of the cells, the
# boundaries of degree 1 (a[] - []a) omitted with :limit, the cells in bar form.
Complexes._basis_title(::CellBasis) = "Cells"  # COV_EXCL_LINE (a constant, folded by the compiler)
Complexes._omits_first_boundaries(::CellBasis) = true  # COV_EXCL_LINE (a constant, folded by the compiler)
Complexes._write_plain_generator(io::IO, w::_CellGenerators, B::CellBasis, i::Integer) =
    _write_cell(io, w.cells, B.tree, first(cell_range(B.tree, B.degree)) + i - 1)

const _LCMMatrix = LabelledSparseMatrix{NTuple{2, Int32}, Int}

"""
    LCMResolution <: GarsideHomology.Complexes.AbstractChainComplex

The lcm resolution of a germ for an ordering, as returned by [`lcm_resolution`](@ref). It is a chain complex of free bimodules over `MonoidAlgebra(ZZ, G, :both)`, which keeps the tree of its cells ([`cell_tree`](@ref)), hence the germ and the ordering.

It implements the interface of `AbstractChainComplex` (`range`, `R[n]`, [`differential`](@ref), `base_ring`) through the `ChainComplex` that it wraps, whose bases are [`CellBasis`](@ref)es. The extension of scalars [`extend_scalars`](@ref) gives a `ChainComplex` over ℤ on the same bases, which still know their cells. The functions that need more than a chain complex, such as the writing of snapshots ([`write_snapshot`](@ref)), have methods for this type.

The cells, the elements and the differential of `R` are obtained as follows:
- its cells, as [`Cell`](@ref)s: `basis(R[n])`, `basis(R)` (the [`GradedBasis`](@ref)), `Cell(R, "2|12")` (from the bar form) and `Cell(R, G(2), G(1, 2))` (from the letters);
- the elements of its modules ([`BasedFreeModuleElement`](@ref)): `R(x)`, the generator `[x]` of a cell `x`, which is also written `R("2|12")` or `R(G(2), G(1, 2))`, and `R("2×1[2]1 - [1]")`, read from its written form;
- its differential: `d = differential(R)`, and `d(x)` for a cell or an element `x` ([`Differential`](@ref)).

At the REPL, the detailed display of `R` gives the cells of each degree and the boundaries `d(x)`, within the screen: it shows the beginning and the end of the list, and omits the boundaries of degree 1, `a[] - []a`, but for a reminder. `show(stdout, "text/plain", R)`, or the same into a file, writes all of them.
"""
struct LCMResolution{G <: AbstractGerm, O <: LCMOrdering, R <: MonoidAlgebra, M <: AbstractModule} <:
       AbstractChainComplex{R, _LCMMatrix}
    tree::CellTree{G, O}
    complex::ChainComplex{R, M, _LCMMatrix}
    function LCMResolution(tree::CellTree{G, O}, complex::ChainComplex{R, M, _LCMMatrix}) where {G, O, R, M}
        A = base_ring(complex)
        A.side === :both && A.monoid === tree.germ ||
            throw(ArgumentError("LCMResolution: the complex is not over the bimodules of the germ of the tree"))
        for n in range(complex)
            b = basis(complex[n])
            b isa CellBasis && b.tree === tree && b.degree == n ||
                throw(ArgumentError("LCMResolution: the basis of degree $n is not that of the cells of the tree"))
        end
        return new{G, O, R, M}(tree, complex)
    end
end

Base.range(R::LCMResolution) = range(R.complex)
ncells(R::LCMResolution) = ncells(R.tree)
ncells(R::LCMResolution, n::Integer) = ncells(R.tree, n)
Base.getindex(R::LCMResolution, n::Integer) = R.complex[n]
Complexes.differential(R::LCMResolution, n::Integer) = differential(R.complex, n)
AbstractAlgebra.base_ring(R::LCMResolution) = base_ring(R.complex)

function Base.show(io::IO, R::LCMResolution)
    T = R.tree
    terms = sum((nnz(differential(R, n)) for n in first(range(R)):(last(range(R)) - 1)); init = 0)
    print(io, "LCMResolution: ", ncells(T), " cells, by degree ", _degree_counts(T), ", ", terms,
          " terms, germ ", nameof(typeof(T.germ)), " with ", natoms(T.germ), " atoms, ")
    show(io, T.ordering)
end

"""
    lcm_resolution(G::AbstractGerm, ordering::LCMOrdering = LCMOrdering())
    lcm_resolution(P::LCMProjectionData)

The lcm resolution of the germ `G` for the ordering `ordering`, as an [`LCMResolution`](@ref), a chain complex over `MonoidAlgebra(ZZ, G, :both)`. It is the free resolution of the monoid algebra `ℤ[M]` by `ℤ[M]`-bimodules obtained from the Garside complex by the matching of the cells (see [`LCMOrdering`](@ref), [`CellTree`](@ref) and [`LCMProjectionData`](@ref)).

The basis of its module of degree `n` is made of the cells of degree `n` ([`CellBasis`](@ref)). The matrix `differential(R, n)`, from degree `n + 1` to degree `n`, is a `LabelledSparseMatrix` whose column `j` is the boundary of the `j`-th cell of degree `n + 1`, with the labels `(σ, τ)` for the terms `c·σ[y]τ`.

Its tensor products by modules of coefficients (`homology(R, M)`, `homology(R, ℤ)`, `ℤ ⊗ R ⊗ ℤ`, [`extend_scalars`](@ref)) compute the homology of the monoid. Since the monoid is a Garside monoid, this homology is also that of its enveloping group, with the same coefficients; this does not hold for a general monoid (see the page Background and conventions of the User guide).

The main uses are:
- `lcm_resolution(G)`, the reduced resolution ([`ComplementDeglex`](@ref));
- `lcm_resolution(G, LCMOrdering(TrivialPreorder(); root_order = p))`, the order resolution of Dehornoy and Lafont for the order `p` of the atoms;
- with Chevie loaded, `lcm_resolution(M, ordering)` for a `BraidMonoid` or a `DualBraidMonoid` `M` of Chevie, which is the resolution of its germ, `lcm_resolution(germ(M), ordering)` ([`germ`](@ref)).

The method with a germ drops the memo of the projection at the end. The method with data `P` reads the boundaries from these data, built by the caller, who keeps the memo.
"""
lcm_resolution(G::AbstractGerm, ordering::LCMOrdering = LCMOrdering()) =
    lcm_resolution(LCMProjectionData(CellTree(G, ordering)))

function lcm_resolution(P::LCMProjectionData)
    T = P.tree
    A = MonoidAlgebra(Complexes.AbstractAlgebra.ZZ, T.germ, :both)
    top = _top_degree(T)
    modules = [BasedFreeModule(A, CellBasis(T, n)) for n in 0:top]
    maps = _LCMMatrix[_differential_matrix(P, n) for n in 0:(top - 1)]
    return LCMResolution(T, ChainComplex(modules, maps; lo = 0))
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
    cell_tree(R::LCMResolution) -> CellTree
    cell_tree(C::AbstractChainComplex) -> CellTree

The tree of the cells of an lcm resolution, or of one of its extensions of scalars, read from the bases of its modules. An `ArgumentError` is thrown if the bases of `C` are not cells.
"""
cell_tree(R::LCMResolution) = R.tree

"""
    (C::AbstractChainComplex)(x::Cell)

The generator `[x]` of the module of the degree of the cell `x`, for a complex whose bases are cells: `P(x)` is `P[n](x)` for a cell `x` of degree `n`.
"""
function (C::AbstractChainComplex)(x::Cell)
    return C[cell_degree(x)](x)
end

"""
    (C::AbstractChainComplex)(letters::GermElement...)

The generator `[x]` of the cell `x` whose letters are `letters`, for a complex whose bases are cells (an lcm resolution, or one of its extensions of scalars). For instance, `P(G(2), G(1, 2))` is the generator `[2|12]`, and `P()` is the generator `[]` of the cell of degree 0. An `ArgumentError` is thrown if the letters are not those of a cell. The cell itself is given by [`Cell`](@ref)`(P, G(2), G(1, 2))`.
"""
(C::AbstractChainComplex)(letters::GermElement...) = C(Cell(C, letters...))

"""
    Cell(C::AbstractChainComplex, text::AbstractString)
    Cell(C::AbstractChainComplex, letters::GermElement...)

The cell of a complex whose bases are cells (an lcm resolution, or one of its extensions of scalars), found from its bar form `text`, read with the settings of the session, or from its letters. For instance, `Cell(P, "2|12")` and `Cell(P, G(2), G(1, 2))` are the cell `2|12`, and `Cell(P, "e")` and `Cell(P)` are the cell of degree 0. An `ArgumentError` is thrown if the text or the letters do not give a cell of the complex.

Calling the complex instead gives the generator of the cell: `P("2|12")` and `P(G(2), G(1, 2))` are the generator `[2|12]`.
"""
Cell(C::AbstractChainComplex, text::AbstractString) = _read_cell(cell_tree(C), notation(), text)

function Cell(C::AbstractChainComplex, letters::GermElement...)
    T = cell_tree(C)
    for s in letters
        s.germ === T.germ || throw(ArgumentError("the letters of a cell are elements of the germ of the resolution"))
    end
    return _cell_of_letters(T, (s.index for s in letters), join(letters, '|'))
end

function cell_tree(C::AbstractChainComplex)
    b = basis(C[first(range(C))])
    b isa CellBasis || throw(ArgumentError("cell_tree: the complex is not an lcm resolution (its bases are not CellBasis)"))
    return b.tree
end
