# The tree of the cells: the cells of the lcm resolution of a germ for an ordering, as a
# prefix tree in integers (J1.d, task d2).

"""
    CellTree(G::AbstractGerm, ordering::LCMOrdering = LCMOrdering())

The cells of the lcm resolution of the germ `G` for the ordering `ordering`, as a tree whose
nodes are the cells. The node of a cell `x = [τ₁|…|τₙ]` has for parent the cell
`[τ₁|…|τₙ₋₁]`, for letter its last letter `τₙ`, for atom the atom `a` that created it
(`τₙ = ω \\ (ω ∨ a)`, `ω` being the weight of the parent), for weight `τ₁⋯τₙ` and for
degree `n`. The node `1` is the cell of degree 0, of weight the identity.

The tree is built breadth first, by the rule of [`LCMOrdering`](@ref): at a node of
weight `ω` whose list `I` of admissible atoms is in the inherited order (the order of the
root at the root), `I` is sorted by the preorder at `ω`, ties keeping the inherited
order; an atom `j` is kept if no atom placed before it divides `ω ∨ aⱼ`; the atoms kept
give the children, in that order, and the list of a child is formed of the atoms of its
elder siblings. No list is copied: it is read in the tree itself.

The cells are numbered `1:ncells(T)` in the order of the construction: by degree, then by
parent, then in the order given by the ordering. The children of a node are contiguous, and the
cells of degree `n` form [`cell_range`](@ref)`(T, n)`, which fixes the order of the basis
of the free module of degree `n`.

Accessors (public, not exported): [`ncells`](@ref), [`cell_range`](@ref),
[`cell_weight`](@ref), [`cell_letter`](@ref), [`cell_atom`](@ref),
[`cell_degree`](@ref), [`cell_children`](@ref), [`admissible_atoms`](@ref),
[`bar_form`](@ref), [`atom_sequence`](@ref).
"""
struct CellTree{G <: AbstractGerm, O <: LCMOrdering}
    germ::G
    ordering::O
    root_order::Vector{Int32}    # the list of the root
    parent::Vector{Int32}        # parent[x], 0 for the root
    atom::Vector{Int32}          # atom[x], 0 for the root
    letter::Vector{Int32}        # letter[x], a simple; 0 for the root
    weight::Vector{Int32}        # weight[x], a simple
    degree::Vector{UInt8}        # degree[x]
    child_start::Vector{Int32}   # children of x: child_start[x]:(child_start[x + 1] - 1)
    degree_start::Vector{Int32}  # degree n: degree_start[n + 1]:(degree_start[n + 2] - 1)
end

# A candidate of a node: an atom of its list, its lcm with the weight, the complement.
const _Candidate = NamedTuple{(:atom, :lcm, :complement), NTuple{3, Int32}}

# The candidate in the form of a result of lcm_atom, as precedes expects it.
_lcm_pair(c::_Candidate) = (lcm = c.lcm, complement = c.complement)

function CellTree(G::AbstractGerm, ordering::LCMOrdering = LCMOrdering())
    root = _root_order(ordering, G)
    parent = Int32[0]
    atom = Int32[0]
    letter = Int32[0]
    weight = Int32[one(G)]
    degree = UInt8[0]
    child_start = Int32[]
    buffer = Vector{_Candidate}(undef, natoms(G))
    B = typeof(left_descents(G, one(G)))
    x = 1
    while x <= length(parent)
        length(parent) < typemax(Int32) || throw(ArgumentError("CellTree: more than $(typemax(Int32)) cells"))
        push!(child_start, Int32(length(parent) + 1))
        # The list of x: the order of the root, or the atoms of the elder siblings of x.
        lo = x == 1 ? 1 : Int(child_start[parent[x]])
        list = x == 1 ? view(root, 1:length(root)) : view(atom, lo:(x - 1))
        _expand!(parent, atom, letter, weight, degree, buffer, G, ordering.preorder, B, x, list)
        x += 1
    end
    push!(child_start, Int32(length(parent) + 1))
    return CellTree(G, ordering, root, parent, atom, letter, weight, degree, child_start,
                    _degree_starts(degree))
end

# The children of the node x, appended to the arrays of the tree. A function of its own,
# so that the weight ω of x is an argument and not a variable of the loop of CellTree.
function _expand!(parent, atom, letter, weight, degree, buffer, G, P, ::Type{B}, x::Int, list) where {B}
    ω = weight[x]
    # 1. The candidates, read from the list before anything is appended to the tree.
    k = 0
    for j in list
        r = lcm_atom(G, ω, j)
        iszero(r.complement) &&
            throw(ArgumentError("CellTree: the atom $j of the list of the cell $x divides its weight $ω"))
        k += 1
        buffer[k] = (atom = Int32(j), lcm = Int32(r.lcm), complement = Int32(r.complement))
    end
    # 2. Stable sort by the preorder at ω: ties keep the inherited order.
    _sort_candidates!(buffer, k, P, G, ω)
    # 3. An atom is kept if no atom placed before it divides its lcm.
    d = degree[x]
    d < typemax(UInt8) || throw(ArgumentError("CellTree: degree above $(typemax(UInt8))"))
    seen = B()
    for i in 1:k
        c = buffer[i]
        if isdisjoint(left_descents(G, c.lcm), seen)
            push!(parent, Int32(x))
            push!(atom, c.atom)
            push!(letter, c.complement)
            push!(weight, c.lcm)
            push!(degree, d + 0x01)
        end
        seen = union(seen, B((c.atom,)))
    end
    return nothing
end

# Insertion sort of buffer[1:k], stable: an element moves before another only if it
# strictly precedes it. The lists have at most N elements.
function _sort_candidates!(buffer, k, P::AbstractLCMPreorder, G, ω)
    for i in 2:k
        c = buffer[i]
        x = _lcm_pair(c)
        j = i - 1
        while j >= 1 && precedes(P, G, ω, x, _lcm_pair(buffer[j]))
            buffer[j + 1] = buffer[j]
            j -= 1
        end
        buffer[j + 1] = c
    end
    return buffer
end

# All lcms tied: the inherited order stays.
_sort_candidates!(buffer, k, ::TrivialPreorder, G, ω) = buffer

# The first cell of each degree, the degrees of the breadth-first order being
# nondecreasing; one more entry, ncells + 1, at the end.
function _degree_starts(degree::Vector{UInt8})
    top = Int(degree[end])
    starts = Vector{Int32}(undef, top + 2)
    starts[1] = 1
    n = 0
    for x in eachindex(degree)
        Int(degree[x]) >= n || throw(ArgumentError("CellTree: degrees not in breadth-first order"))
        while Int(degree[x]) > n
            n += 1
            starts[n + 1] = x
        end
    end
    starts[top + 2] = length(degree) + 1
    return starts
end

# ------------------------------------------------------------------------------------
# Accessors

"""
    ncells(T::CellTree) -> Int
    ncells(T::CellTree, n) -> Int
    ncells(R::LCMResolution) -> Int
    ncells(R::LCMResolution, n) -> Int

The number of cells of `T`, or of its cells of degree `n` (0 outside the degrees of `T`);
for an lcm resolution `R`, those of its tree ([`cell_tree`](@ref)). The numbers by degree:
`[ncells(R, n) for n in range(R)]`.
"""
ncells(T::CellTree) = length(T.parent)
ncells(T::CellTree, n::Integer) = length(cell_range(T, n))

# The largest degree of a cell of T.
_top_degree(T::CellTree) = length(T.degree_start) - 2

# The numbers of cells by degree, from 0.
_degree_counts(T::CellTree) = [ncells(T, n) for n in 0:_top_degree(T)]

"""
    cell_range(T::CellTree, n) -> UnitRange{Int}

The numbers of the cells of degree `n` (empty outside the degrees of `T`), in the order
of the basis of the free module of degree `n`.
"""
cell_range(T::CellTree, n::Integer) =
    0 <= n <= _top_degree(T) ? (Int(T.degree_start[n + 1]):Int(T.degree_start[n + 2] - 1)) : (1:0)

"""
    cell_weight(T::CellTree, x) -> Int32

The weight of the cell `x`, the product of its letters, a simple of the germ.
"""
cell_weight(T::CellTree, x::Integer) = T.weight[x]

"""
    cell_letter(T::CellTree, x) -> Int32

The last letter of the cell `x` (a simple), `0` for the cell of degree 0.
"""
cell_letter(T::CellTree, x::Integer) = T.letter[x]

"""
    cell_atom(T::CellTree, x) -> Int32

The atom that created the cell `x` from its parent, `0` for the cell of degree 0.
"""
cell_atom(T::CellTree, x::Integer) = T.atom[x]

"""
    cell_degree(T::CellTree, x) -> Int

The degree of the cell `x`.
"""
cell_degree(T::CellTree, x::Integer) = Int(T.degree[x])

"""
    cell_children(T::CellTree, x) -> UnitRange{Int}

The children of the cell `x`, contiguous, in the order given by the ordering.
"""
cell_children(T::CellTree, x::Integer) = Int(T.child_start[x]):Int(T.child_start[x + 1] - 1)

"""
    admissible_atoms(T::CellTree, x)

The list of admissible atoms of the cell `x`, in the inherited order, as a view of the
tree: the order of the root for the cell of degree 0, the atoms of the elder siblings of
`x` otherwise.
"""
admissible_atoms(T::CellTree, x::Integer) =
    x == 1 ? view(T.root_order, 1:length(T.root_order)) : view(T.atom, Int(T.child_start[T.parent[x]]):(Int(x) - 1))

"""
    bar_form(T::CellTree, x) -> Vector{Int32}

The letters of the cell `x`, from the first: the simples `[τ₁, …, τₙ]`.
"""
bar_form(T::CellTree, x::Integer) = _path(T.letter, T, x)

"""
    atom_sequence(T::CellTree, x) -> Vector{Int32}

The atoms that created the cell `x` and its ancestors, from the first: the name of the
cell in the snapshots (`7.3.12`).
"""
atom_sequence(T::CellTree, x::Integer) = _path(T.atom, T, x)

function _path(field::Vector{Int32}, T::CellTree, x::Integer)
    n = cell_degree(T, x)
    w = Vector{Int32}(undef, n)
    y = Int(x)
    for i in n:-1:1
        w[i] = field[y]
        y = Int(T.parent[y])
    end
    return w
end

function Base.show(io::IO, T::CellTree)
    print(io, "CellTree: ", ncells(T), " cells, by degree ", _degree_counts(T), ", ")
    show(io, T.ordering)
end
