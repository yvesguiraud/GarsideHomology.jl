# The tree of the chains: the chains of the Garside complex visited by the projection,
# as a prefix tree whose first nodes are the cells (J1.d, task d3).

"""
    ChainTree(T::CellTree)

The prefix tree of the chains `[x₁|…|xₙ]` of the Garside complex visited by the
projection of a resolution by lcms: the node of a chain has for parent the chain
`[x₁|…|xₙ₋₁]`, for letter `xₙ` (a simple other than the identity) and for degree `n`.
Its first nodes are the cells of `T`, with the same numbers, so that a chain is a cell if
and only if its number is at most `ncells(T)`; the other chains are added as they are
met. The child of a node by a letter is found in a dictionary, whose key packs the
number of the node and the letter in a `UInt64`.
"""
struct ChainTree
    parent::Vector{Int32}
    letter::Vector{Int32}
    degree::Vector{UInt8}
    index::Dict{UInt64, Int32}
    ncells::Int
end

_chain_key(p::Integer, s::Integer) = (UInt64(UInt32(p)) << 32) | UInt64(UInt32(s))

function ChainTree(T::CellTree)
    n = ncells(T)
    index = Dict{UInt64, Int32}()
    sizehint!(index, n)
    for x in 2:n
        key = _chain_key(T.parent[x], T.letter[x])
        haskey(index, key) && throw(ArgumentError("ChainTree: two children of the cell $(T.parent[x]) with the letter $(T.letter[x])"))
        index[key] = Int32(x)
    end
    return ChainTree(copy(T.parent), copy(T.letter), copy(T.degree), index, n)
end

# The number of chains of the tree.
_nchains(C::ChainTree) = length(C.parent)

# The child of the node p by the letter s, 0 if it is not in the tree.
_child(C::ChainTree, p::Integer, s::Integer) = get(C.index, _chain_key(p, s), Int32(0))

# The child of the node p by the letter s, added to the tree if needed.
function _child!(C::ChainTree, p::Integer, s::Integer)
    iszero(s) && throw(ArgumentError("ChainTree: the identity is not a letter of a chain"))
    key = _chain_key(p, s)
    x = get(C.index, key, Int32(0))
    iszero(x) || return x
    d = C.degree[p]
    d < typemax(UInt8) || throw(ArgumentError("ChainTree: degree above $(typemax(UInt8))"))
    push!(C.parent, Int32(p))
    push!(C.letter, Int32(s))
    push!(C.degree, d + 0x01)
    x = Int32(length(C.parent))
    C.index[key] = x
    return x
end

# The letters of the chain x, from the first.
function _letters(C::ChainTree, x::Integer)
    n = Int(C.degree[x])
    w = Vector{Int32}(undef, n)
    y = Int(x)
    for i in n:-1:1
        w[i] = C.letter[y]
        y = Int(C.parent[y])
    end
    return w
end

# The prefixes of the chain x: a[k + 1] is the node of its prefix of length k, k = 0:n − 1.
function _prefixes(C::ChainTree, x::Integer)
    n = Int(C.degree[x])
    a = Vector{Int32}(undef, max(n, 1))
    y = Int32(x)
    for k in n:-1:1
        y = C.parent[y]
        a[k] = y
    end
    n == 0 && (a[1] = Int32(1))
    return a
end

# The node of the chain obtained from the node x by appending the letters letters[r].
function _extend!(C::ChainTree, x::Integer, letters, r)
    y = Int32(x)
    for k in r
        y = _child!(C, y, letters[k])
    end
    return y
end

function Base.show(io::IO, C::ChainTree)
    print(io, "ChainTree: ", _nchains(C), " chains, of which ", C.ncells, " cells")
end
