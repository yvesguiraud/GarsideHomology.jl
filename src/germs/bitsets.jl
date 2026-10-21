# Fixed-width bit sets, used for the sets of atoms of a germ (left descents).

# Marker of the inner constructor, which takes the words themselves.
struct _Words end

"""
    StaticBitSet{K}

Set of integers of `1:64K`, stored as `K` words of 64 bits: the integer `i` is the bit
`i - 1`. The value is immutable and isbits, so that sets are stored inline in arrays
and returned without allocation; `K` is a constant of the type, fixed for instance by a
germ from its number of atoms (`K = cld(N, 64)`).

`StaticBitSet{K}()` is the empty set and `StaticBitSet{K}(itr)` the set of the integers
of `itr`, each in `1:64K`. Fast methods: `in`, iteration in increasing order, `length`,
`isempty`, `first`, `issubset`, `union`, `intersect`, `setdiff`, `symdiff`,
`isdisjoint`, `==`. As an `AbstractSet{Int}`, a `StaticBitSet` also has the generic
methods of sets (`collect`, comparison and union with other sets, `hash`).
"""
struct StaticBitSet{K} <: AbstractSet{Int}
    words::NTuple{K, UInt64}
    function StaticBitSet{K}(::_Words, words::NTuple{K, UInt64}) where {K}
        K isa Int && K >= 1 || throw(ArgumentError("StaticBitSet{$K}: K must be an Int ≥ 1"))
        return new{K}(words)
    end
end

_bitset(words::NTuple{K, UInt64}) where {K} = StaticBitSet{K}(_Words(), words)

StaticBitSet{K}() where {K} = _bitset(ntuple(_ -> zero(UInt64), Val(K)))

function StaticBitSet{K}(itr) where {K}
    words = ntuple(_ -> zero(UInt64), Val(K))
    for i in itr
        i isa Integer || throw(ArgumentError("StaticBitSet{$K}: element $i is not an integer"))
        j = Int(i) - 1
        0 <= j < 64K || throw(ArgumentError("StaticBitSet{$K}: element $i out of 1:$(64K)"))
        k = (j >> 6) + 1
        words = Base.setindex(words, words[k] | (one(UInt64) << (j & 63)), k)
    end
    return _bitset(words)
end

@inline function Base.in(i::Integer, s::StaticBitSet{K}) where {K}
    j = Int(i) - 1
    0 <= j < 64K || return false
    @inbounds w = s.words[(j >> 6) + 1]
    return (w >> (j & 63)) & one(UInt64) != 0
end

@inline function Base.iterate(s::StaticBitSet{K}, state::Tuple{Int, UInt64} = (1, s.words[1])) where {K}
    k, w = state
    while w == 0
        k += 1
        k > K && return nothing
        @inbounds w = s.words[k]
    end
    return (64 * (k - 1) + trailing_zeros(w) + 1, (k, w & (w - one(UInt64))))
end

Base.length(s::StaticBitSet) = sum(count_ones, s.words)
Base.isempty(s::StaticBitSet) = all(iszero, s.words)

function Base.first(s::StaticBitSet)
    r = iterate(s)
    r === nothing && throw(ArgumentError("collection must be non-empty"))
    return r[1]
end

Base.:(==)(s::StaticBitSet{K}, t::StaticBitSet{K}) where {K} = s.words == t.words
Base.issubset(s::StaticBitSet{K}, t::StaticBitSet{K}) where {K} =
    all(iszero, map((a, b) -> a & ~b, s.words, t.words))
Base.isdisjoint(s::StaticBitSet{K}, t::StaticBitSet{K}) where {K} =
    all(iszero, map(&, s.words, t.words))

Base.union(s::StaticBitSet{K}, ts::StaticBitSet{K}...) where {K} =
    _bitset(foldl((a, t) -> map(|, a, t.words), ts; init = s.words))
Base.intersect(s::StaticBitSet{K}, ts::StaticBitSet{K}...) where {K} =
    _bitset(foldl((a, t) -> map(&, a, t.words), ts; init = s.words))
Base.setdiff(s::StaticBitSet{K}, t::StaticBitSet{K}) where {K} =
    _bitset(map((a, b) -> a & ~b, s.words, t.words))
Base.symdiff(s::StaticBitSet{K}, t::StaticBitSet{K}) where {K} =
    _bitset(map(xor, s.words, t.words))

Base.copy(s::StaticBitSet) = s

function Base.show(io::IO, s::StaticBitSet{K}) where {K}
    print(io, "StaticBitSet{", K, "}(")
    show(io, collect(s))
    print(io, ")")
end
