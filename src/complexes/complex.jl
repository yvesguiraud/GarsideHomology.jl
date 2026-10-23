# Chain complexes of free modules.

"""
    AbstractChainComplex{R, D}

Supertype of the chain complexes `C_lo ← C_{lo+1} ← … ← C_hi` of this layer, over a ring
described by an object of type `R` (an AbstractAlgebra ring or a [`MonoidAlgebra`](@ref)),
with maps of type `D` (`SparseMatrixCSC` over ℤ, [`LabelledSparseMatrix`](@ref) over a
monoid algebra). Convention of Gauss: `d_n : C_{n+1} → C_n` has the index of its target.

Interface, required: `range(C)` (the degrees `lo:hi`), `C[n]` (the module `C_n`, a
`BoundsError` outside the degrees), [`differential`](@ref)`(C, n)` (the map `d_n`, the zero
map outside `lo:(hi - 1)`) and `base_ring(C)`. Derived, for every subtype: `firstindex`,
`lastindex`, [`specialise`](@ref), [`is_square_zero`](@ref), [`homology`](@ref),
[`collapsible_pairs`](@ref). [`ChainComplex`](@ref) is the generic concrete type; a type
that carries more structure (the lcm resolutions of Gauss, for instance) is a subtype of
its own, usually wrapping a `ChainComplex`.
"""
abstract type AbstractChainComplex{R <: _RingDescriptor, D} end

Base.firstindex(C::AbstractChainComplex) = first(range(C))
Base.lastindex(C::AbstractChainComplex) = last(range(C))

"""
    ChainComplex(modules, maps; lo = 0, check = true)
    ChainComplex(maps::AbstractVector{<:SparseMatrixCSC}; lo = 0)

Chain complex `C_lo ← C_{lo+1} ← … ← C_hi` of free modules over a common ring, where
`maps[k]` is `d_n : C_{n+1} → C_n` for `n = lo + k - 1` (convention of Gauss: `d_n` has
target `C_n`, whereas Hecke and most texts index `d_n : C_n → C_{n-1}` by the source).
The maps are [`LabelledSparseMatrix`](@ref)s over a [`MonoidAlgebra`](@ref), or
`SparseMatrixCSC`s over ℤ. The second form takes integer matrices only and builds the
modules `BasedFreeModule(ZZ, r)` of the ranks they imply.

Checked: one map fewer than modules, a common ring, `size(d_n) == (rank(C_n), rank(C_{n+1}))`.
Not checked: `d_n ∘ d_{n+1} = 0` (see [`is_square_zero`](@ref)).

Accessors: `range(C)` (degrees `lo:hi`), `C[n]` (the module `C_n`),
[`differential`](@ref)`(C, n)` (the map `d_n`), `base_ring(C)`.
"""
struct ChainComplex{R <: _RingDescriptor, M <: AbstractModule, D} <: AbstractChainComplex{R, D}
    ring::R
    modules::Vector{M}
    maps::Vector{D}
    lo::Int
    function ChainComplex{R, M, D}(ring, modules::Vector{M}, maps::Vector{D}, lo::Integer;
                                   check::Bool = true) where {R, M, D}
        C = new{R, M, D}(ring, modules, maps, Int(lo))
        check && _check_complex(C)
        return C
    end
end

function ChainComplex(modules::AbstractVector{M}, maps::AbstractVector{D};
                      lo::Integer = 0, check::Bool = true) where {M <: AbstractModule, D}
    isempty(modules) && throw(ArgumentError("a chain complex needs at least one module"))
    ring = base_ring(first(modules))
    return ChainComplex{typeof(ring), M, D}(ring, collect(modules), collect(maps), lo; check)
end

function ChainComplex(maps::AbstractVector{<:SparseMatrixCSC{<:Integer}}; lo::Integer = 0)
    isempty(maps) &&
        throw(ArgumentError("give at least one matrix, or use ChainComplex(modules, maps)"))
    ranks = [size(first(maps), 1); [size(A, 2) for A in maps]]
    modules = [BasedFreeModule(AbstractAlgebra.ZZ, r) for r in ranks]
    return ChainComplex(modules, collect(maps); lo)
end

function _check_complex(C::ChainComplex)
    nm, nd = length(C.modules), length(C.maps)
    nm >= 1 || throw(ArgumentError("a chain complex needs at least one module"))
    nd == nm - 1 || throw(ArgumentError("$nm modules need $(nm - 1) maps, got $nd"))
    for (k, M) in enumerate(C.modules)
        base_ring(M) == C.ring ||
            throw(ArgumentError("the module in degree $(C.lo + k - 1) is over another ring"))
    end
    for k in 1:nd
        n = C.lo + k - 1
        expected = (rank(C.modules[k]), rank(C.modules[k + 1]))
        size(C.maps[k]) == expected ||
            throw(ArgumentError("d_$n has size $(size(C.maps[k], 1))×$(size(C.maps[k], 2)), " *
                                "expected rank(C_$n) × rank(C_$(n + 1)) = $(expected[1])×$(expected[2])"))
    end
    return C
end

Base.range(C::ChainComplex) = C.lo:(C.lo + length(C.modules) - 1)
AbstractAlgebra.base_ring(C::ChainComplex) = C.ring

function Base.getindex(C::ChainComplex, n::Integer)
    n in range(C) || throw(BoundsError(C, n))
    return C.modules[n - C.lo + 1]
end

_rank(C::AbstractChainComplex, n::Integer) = n in range(C) ? rank(C[n]) : 0

"""
    differential(C::AbstractChainComplex, n)

The map `d_n : C_{n+1} → C_n`, indexed by its target (see [`AbstractChainComplex`](@ref)), of
size `rank(C_n) × rank(C_{n+1})`. Outside the stored
maps it is the zero map, of size `0 × rank(C_lo)` for `n = lo - 1` and
`rank(C_hi) × 0` for `n = hi`, as the formula of homology expects.
"""
function differential(C::ChainComplex{R, M, D}, n::Integer) where {R, M, D}
    k = n - C.lo + 1
    1 <= k <= length(C.maps) && return C.maps[k]
    return _zero_map(D, _rank(C, n), _rank(C, n + 1))
end

_zero_map(::Type{SparseMatrixCSC{Tv, Ti}}, m, n) where {Tv, Ti} = spzeros(Tv, Ti, m, n)
_zero_map(::Type{LabelledSparseMatrix{L, Tv}}, m, n) where {L, Tv} = LabelledSparseMatrix{L, Tv}(m, n)

"""
    specialise(C::AbstractChainComplex)
    specialise(C::AbstractChainComplex, χ)

The [`ChainComplex`](@ref) over the ring of coefficients of the monoid algebra of `C`, with
the same bases (shared, not copied) and the maps `specialise(d_n)` or `specialise(d_n, χ)`.
For a resolution of Gauss, the typed interface is `homology(R, M)` and
`Gauss.tensor_product`, with a module of coefficients `M`.
"""
function specialise(C::AbstractChainComplex{<:MonoidAlgebra, LabelledSparseMatrix{L, Tv}},
                    χ = nothing) where {L, Tv}
    base = base_ring(base_ring(C))
    degrees = range(C)
    modules = [BasedFreeModule(base, cells(C[n])) for n in degrees]
    maps = SparseMatrixCSC{Tv, Int}[_specialise(differential(C, n), χ) for n in first(degrees):(last(degrees) - 1)]
    return ChainComplex(modules, maps; lo = first(degrees))
end

"""
    is_square_zero(C::AbstractChainComplex)

Whether `d_n ∘ d_{n+1} = 0` for every `n`, for a complex of integer matrices, with
checked arithmetic. Over a monoid algebra this needs products of elements of the monoid:
specialise the complex first, or use the certificates of the resolution layer.
"""
function is_square_zero(C::AbstractChainComplex{R, D}) where {R, D <: SparseMatrixCSC{<:Integer}}
    degrees = range(C)
    for n in first(degrees):(last(degrees) - 2)
        _product_is_zero(differential(C, n), differential(C, n + 1)) || return false
    end
    return true
end

is_square_zero(::AbstractChainComplex{R, D}) where {R, D <: LabelledSparseMatrix} =
    throw(ArgumentError("d² = 0 over a monoid algebra needs products of monoid elements: " *
                        "specialise the complex, or use the certificates of the resolution"))

# Whether A * B = 0, one column of B at a time, with a dense accumulator.
function _product_is_zero(A::SparseMatrixCSC{Ta}, B::SparseMatrixCSC{Tb}) where {Ta, Tb}
    T = promote_type(Ta, Tb)
    acc = zeros(T, size(A, 1))
    touched = Int[]
    Ar, Av, Br, Bv = rowvals(A), nonzeros(A), rowvals(B), nonzeros(B)
    for j in 1:size(B, 2)
        for q in nzrange(B, j)
            l, b = Br[q], T(Bv[q])
            for p in nzrange(A, l)
                i = Ar[p]
                iszero(acc[i]) && push!(touched, i)
                acc[i] = _checked_add(acc[i], _checked_mul(T(Av[p]), b))
            end
        end
        for i in touched
            iszero(acc[i]) || return false
        end
        empty!(touched)
    end
    return true
end

function Base.show(io::IO, C::ChainComplex)
    print(io, "ChainComplex over ", C.ring, " in degrees ", range(C), ", ranks ",
          join((rank(N) for N in C.modules), ", "))
end
