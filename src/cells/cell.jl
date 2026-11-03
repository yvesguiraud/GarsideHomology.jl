# The cells for the user (J1.h): a cell of the tree, written in bar form x₁|⋯|xₙ, read as the
# sequence of its letters, and found again from its letters or its written form. The tree
# keeps its numbers (Int32); only what is handed to the user is wrapped.

"""
    Cell

A cell of an lcm resolution, for the user: the tree of the cells ([`CellTree`](@ref)) and
the number of the cell in it. It is written in bar form, `x₁|⋯|xₙ`, its letters being
written as words in the atoms ([`GermElement`](@ref)), and `e` for the cell of degree 0,
following the settings of the display ([`Notation`](@ref)): `2|12` in B⁺(A₂).

The elements of the bases of the modules of a resolution
([`basis`](@ref)`(P[n])`, a [`CellBasis`](@ref)), and of its extensions of scalars. A cell
is also given by its bar form, `P("2|12")`, or by its letters, `P(G(2), G(1, 2))`.

A cell reads as the sequence of its letters: `c[i]` is the `i`-th letter (a
`GermElement`), `c[end]` the last one, `length(c)` the degree, and `collect(c)` and the
iteration give the letters. Also [`cell_weight`](@ref)`(c)` (the product of the letters),
[`cell_degree`](@ref)`(c)` and [`atom_sequence`](@ref)`(c)` (the atoms that created the
cell, its name in the snapshots). Two cells are equal if they are the same cell of the same
tree; `isless` is the order of the construction. A `Cell` is a scalar for broadcasting.
"""
struct Cell{T <: CellTree}
    tree::T
    index::Int32
end

Base.length(c::Cell) = cell_degree(c.tree, c.index)
Base.firstindex(c::Cell) = 1  # COV_EXCL_LINE (a constant, folded by the compiler)
Base.lastindex(c::Cell) = length(c)
Base.eltype(::Type{<:Cell{<:CellTree{G}}}) where {G} = GermElement{G}

function Base.getindex(c::Cell, i::Integer)
    n = length(c)
    1 <= i <= n || throw(BoundsError(c, i))
    return GermElement(c.tree.germ, bar_form(c.tree, c.index)[i]; check = false)
end

function Base.iterate(c::Cell, state = (bar_form(c.tree, c.index), 1))
    letters, i = state
    i > length(letters) && return nothing
    return GermElement(c.tree.germ, letters[i]; check = false), (letters, i + 1)
end

Base.broadcastable(c::Cell) = Ref(c)

"""
    cell_weight(c::Cell) -> GermElement

The weight of the cell `c`, the product of its letters, as an element of the germ.
"""
cell_weight(c::Cell) = GermElement(c.tree.germ, cell_weight(c.tree, c.index); check = false)

"""
    cell_degree(c::Cell) -> Int

The degree of the cell `c`, the number of its letters.
"""
cell_degree(c::Cell) = cell_degree(c.tree, c.index)

"""
    atom_sequence(c::Cell) -> Vector{GermElement}

The atoms that created the cell `c` and its ancestors, from the first, as elements of the
germ: the name of the cell in the snapshots (`7.3.12`).
"""
atom_sequence(c::Cell) = [GermElement(c.tree.germ, a; check = false) for a in atom_sequence(c.tree, c.index)]

Base.:(==)(c::Cell, d::Cell) = c.tree === d.tree && c.index == d.index
Base.hash(c::Cell, h::UInt) = hash(c.index, hash(objectid(c.tree), hash(:Cell, h)))

function Base.isless(c::Cell, d::Cell)
    c.tree === d.tree || throw(ArgumentError("isless: the two cells belong to different trees"))
    return c.index < d.index
end

# ------------------------------------------------------------------------------------
# Writing

# What the writing of the cells of a tree needs, computed once for a display.
struct _CellWriter{W <: _SimpleWriter}
    simple::W
    bar_separator::String
    empty_cell::String
end

_CellWriter(T::CellTree, N::Notation) = _CellWriter(_SimpleWriter(T.germ, N), N.bar_separator, N.empty_cell)

function _write_cell(io::IO, w::_CellWriter, T::CellTree, x::Integer)
    cell_degree(T, x) == 0 && return print(io, w.empty_cell)
    for (i, s) in enumerate(bar_form(T, x))
        i > 1 && print(io, w.bar_separator)
        _write_simple(io, w.simple, s)
    end
    return nothing
end

Base.show(io::IO, c::Cell) = _write_cell(io, _CellWriter(c.tree, notation(io)), c.tree, c.index)

# ------------------------------------------------------------------------------------
# Finding a cell from its letters or its bar form

# The cell of T whose letters are `letters` (numbers of simples), found by descending the
# tree: the children of a cell have distinct letters, since an atom whose lcm with the
# weight is that of an elder sibling is not kept. `context` names the request in the message.
function _cell_of_letters(T::CellTree, letters, context::AbstractString)
    x = 1
    for s in letters
        found = 0
        for y in cell_children(T, x)
            if T.letter[y] == s
                found = y
                break
            end
        end
        found == 0 && throw(ArgumentError("$context is not a cell of the resolution"))
        x = found
    end
    return Cell(T, Int32(x))
end

# The cell of T written `text` with the settings N.
function _read_cell(T::CellTree, N::Notation, text::AbstractString)
    t = strip(text)
    t == N.empty_cell && return Cell(T, Int32(1))
    parts = split(t, N.bar_separator)
    letters = Int32[_parse_simple(T.germ, N, p; context = text) for p in parts]
    return _cell_of_letters(T, letters, "\"$text\"")
end
