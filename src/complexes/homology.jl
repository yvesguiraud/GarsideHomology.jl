# Finitely generated graded modules, and the homology of complexes over ℤ.

"""
    FinitelyGeneratedGradedModule{T}(modules; lo = 0)
    FinitelyGeneratedGradedAbelianGroup(modules; lo = 0)

Finitely generated graded module `M = ⊕ₙ Mₙ` over a principal ideal domain, concentrated
in the degrees `lo:(lo + length(modules) - 1)`, where `Mₙ = modules[n - lo + 1]` is a
[`FinitelyGeneratedModule{T}`](@ref). `FinitelyGeneratedGradedAbelianGroup` is the alias
for `T = BigInt`, that is over ℤ; it is what [`homology`](@ref) returns.

`M[n]` is the module of degree `n` (a `BoundsError` outside the degrees, as for
`ChainComplex`), `range(M)` the degrees.
"""
struct FinitelyGeneratedGradedModule{T}
    lo::Int
    modules::Vector{FinitelyGeneratedModule{T}}
    function FinitelyGeneratedGradedModule{T}(modules::AbstractVector; lo::Integer = 0) where {T}
        return new{T}(Int(lo), Vector{FinitelyGeneratedModule{T}}(modules))
    end
end

"""
    FinitelyGeneratedGradedAbelianGroup

Alias of `FinitelyGeneratedGradedModule{BigInt}`: a finitely generated graded abelian group.
"""
const FinitelyGeneratedGradedAbelianGroup = FinitelyGeneratedGradedModule{BigInt}

Base.range(M::FinitelyGeneratedGradedModule) = M.lo:(M.lo + length(M.modules) - 1)
Base.firstindex(M::FinitelyGeneratedGradedModule) = M.lo
Base.lastindex(M::FinitelyGeneratedGradedModule) = M.lo + length(M.modules) - 1

function Base.getindex(M::FinitelyGeneratedGradedModule, n::Integer)
    n in range(M) || throw(BoundsError(M, n))
    return M.modules[n - M.lo + 1]
end

Base.:(==)(M::FinitelyGeneratedGradedModule, N::FinitelyGeneratedGradedModule) =
    M.lo == N.lo && M.modules == N.modules
Base.hash(M::FinitelyGeneratedGradedModule, h::UInt) =
    hash(M.modules, hash(M.lo, hash(:FinitelyGeneratedGradedModule, h)))

function Base.show(io::IO, M::FinitelyGeneratedGradedAbelianGroup)
    print(io, "FinitelyGeneratedGradedAbelianGroup in degrees ", range(M), ": ",
          join((string(H) for H in M.modules), "; "))
end

function Base.show(io::IO, ::MIME"text/plain", M::FinitelyGeneratedGradedAbelianGroup)
    print(io, "FinitelyGeneratedGradedAbelianGroup in degrees ", range(M))
    for n in range(M)
        print(io, "\n  ", n, ": ", M[n])
    end
end

# Summary of the elementary divisors of an integer matrix: its rank, the number of
# divisors equal to 1 and the other divisors, in increasing order (internal).
struct _DivisorSummary
    rank::Int
    ones::Int
    torsion::Vector{BigInt}
end

# Matrices without any nonzero entry are handled here; the others go to the solver,
# which passes small or dense cores to Hecke itself.
function _divisor_summary(A::SparseMatrixCSC{<:Integer}; kwargs...)
    (size(A, 1) == 0 || size(A, 2) == 0 || all(iszero, nonzeros(A))) &&
        return _DivisorSummary(0, 0, BigInt[])
    d = sparse_elementary_divisors(A; kwargs...)
    u = count(isone, d)
    return _DivisorSummary(length(d), u, d[(u + 1):end])
end

"""
    homology(C::ChainComplex)
    homology(C::ChainComplex, n)

Homology of a complex over ℤ: a [`FinitelyGeneratedGradedAbelianGroup`](@ref) over the
degrees of `C`, or the [`FinitelyGeneratedAbelianGroup`](@ref) `H_n` alone. With
`d_n : C_{n+1} → C_n`, the free rank of `H_n` is `rank(C_n) − rank(d_n) − rank(d_{n-1})`
and its torsion is given by the elementary divisors of `d_n` other than 1, computed by
`SparseElementaryDivisors.sparse_elementary_divisors`, to which the keyword arguments are
passed (`value_type`, `nworkers`, `core_format`, `verbose`, ...). `homology(C)` computes
each differential once. A complex over a monoid algebra must be specialised first, e.g.
`homology(specialise(C))`.
"""
function homology(C::ChainComplex{R, M, D}; kwargs...) where
        {R <: AbstractAlgebra.Integers, M, D <: SparseMatrixCSC{<:Integer}}
    degrees = range(C)
    summaries = Dict(n => _divisor_summary(differential(C, n); kwargs...)
                     for n in (first(degrees) - 1):last(degrees))
    groups = [FinitelyGeneratedAbelianGroup(rank(C[n]) - summaries[n].rank - summaries[n - 1].rank,
                                            summaries[n].torsion) for n in degrees]
    return FinitelyGeneratedGradedAbelianGroup(groups; lo = first(degrees))
end

function homology(C::ChainComplex{R, M, D}, n::Integer; kwargs...) where
        {R <: AbstractAlgebra.Integers, M, D <: SparseMatrixCSC{<:Integer}}
    n in range(C) || throw(BoundsError(C, n))
    a = _divisor_summary(differential(C, n); kwargs...)
    b = _divisor_summary(differential(C, n - 1); kwargs...)
    return FinitelyGeneratedAbelianGroup(rank(C[n]) - a.rank - b.rank, a.torsion)
end

homology(::ChainComplex{<:MonoidAlgebra}, args...; kwargs...) =
    throw(ArgumentError("homology needs a complex over ℤ: specialise the complex first, " *
                        "e.g. homology(specialise(C))"))
