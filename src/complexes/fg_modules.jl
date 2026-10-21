# Finitely generated modules over a principal ideal domain, as results of homology.

"""
    FinitelyGeneratedModule{T}(rank, torsion = T[])
    FinitelyGeneratedAbelianGroup(rank, torsion = BigInt[])

Finitely generated module `R^rank ⊕ R/(t₁) ⊕ … ⊕ R/(tₖ)` over a principal ideal domain
`R`, stored by its free rank and its invariant factors `t₁ | t₂ | … | tₖ`, none of them a
unit; `T` is the type of the invariant factors. Only `T = BigInt` is implemented for
now, that is `R = ℤ`, under the alias `FinitelyGeneratedAbelianGroup`; polynomial
factors, for `k[t, t⁻¹]`, will follow.

Over ℤ the constructor accepts any list of positive integers: factors equal to 1 are
dropped and the others are brought to invariant-factor form, so that `ℤ/2ℤ ⊕ ℤ/3ℤ`
becomes `ℤ/6ℤ`. Displayed as in Garnier's tables, torsion first: `(ℤ/3ℤ)² ⊕ ℤ/6ℤ ⊕ ℤ`,
and `0` for the trivial group. `parse(FinitelyGeneratedAbelianGroup, s)` reads this form,
and also `ℤ/2`, `Z` for `ℤ`, `^k` for an exponent and `+` for `⊕`.

`rank(H)` is the free rank, [`torsion`](@ref)`(H)` the invariant factors.
"""
struct FinitelyGeneratedModule{T}
    rank::Int
    torsion::Vector{T}
    function FinitelyGeneratedModule{T}(rank::Integer, torsion::AbstractVector) where {T}
        rank >= 0 || throw(ArgumentError("negative rank $rank"))
        return new{T}(Int(rank), _invariant_factors(T, torsion))
    end
end

FinitelyGeneratedModule{T}(rank::Integer) where {T} = FinitelyGeneratedModule{T}(rank, T[])

"""
    FinitelyGeneratedAbelianGroup

Alias of `FinitelyGeneratedModule{BigInt}`: a finitely generated abelian group.
"""
const FinitelyGeneratedAbelianGroup = FinitelyGeneratedModule{BigInt}

# Invariant factors of diag(t) over ℤ: replacing (tᵢ, tⱼ) by their gcd and lcm for all
# i < j yields a divisibility chain (Smith normal form of a diagonal matrix).
function _invariant_factors(::Type{BigInt}, torsion::AbstractVector)
    t = BigInt[]
    for x in torsion
        x isa Integer || throw(ArgumentError("torsion factors over ℤ must be integers, got $x"))
        x >= 1 || throw(ArgumentError("torsion factors must be positive, got $x"))
        isone(x) || push!(t, BigInt(x))
    end
    k = length(t)
    for i in 1:k, j in (i + 1):k
        a, b = t[i], t[j]
        g = gcd(a, b)
        t[i] = g
        t[j] = div(a, g) * b
    end
    return filter!(!isone, t)
end

_invariant_factors(::Type{T}, torsion::AbstractVector) where {T} =
    throw(ArgumentError("invariant factors of type $T are not supported yet (only BigInt, for ℤ)"))

AbstractAlgebra.rank(H::FinitelyGeneratedModule) = H.rank

"""
    torsion(H::FinitelyGeneratedModule)

The invariant factors `t₁ | t₂ | … | tₖ` of the torsion submodule of `H`, none of them a
unit (all at least 2 over ℤ).
"""
torsion(H::FinitelyGeneratedModule) = H.torsion

Base.iszero(H::FinitelyGeneratedModule) = H.rank == 0 && isempty(H.torsion)
Base.:(==)(G::FinitelyGeneratedModule, H::FinitelyGeneratedModule) =
    G.rank == H.rank && G.torsion == H.torsion
Base.hash(H::FinitelyGeneratedModule, h::UInt) =
    hash(H.torsion, hash(H.rank, hash(:FinitelyGeneratedModule, h)))

const _SUPERSCRIPTS = ('⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹')

_superscript(n::Integer) = join(_SUPERSCRIPTS[d + 1] for d in reverse(digits(n)))

function Base.show(io::IO, H::FinitelyGeneratedAbelianGroup)
    parts = String[]
    t = H.torsion
    i = 1
    while i <= length(t)
        j = i
        while j < length(t) && t[j + 1] == t[i]
            j += 1
        end
        e = j - i + 1
        push!(parts, e == 1 ? "ℤ/$(t[i])ℤ" : "(ℤ/$(t[i])ℤ)$(_superscript(e))")
        i = j + 1
    end
    H.rank > 0 && push!(parts, H.rank == 1 ? "ℤ" : "ℤ$(_superscript(H.rank))")
    print(io, isempty(parts) ? "0" : join(parts, " ⊕ "))
end

function Base.parse(::Type{FinitelyGeneratedAbelianGroup}, s::AbstractString)
    str = replace(strip(s), "ℤ" => "Z", " " => "", "⊕" => "+")
    # Runs of superscript digits become ^digits.
    buf = IOBuffer()
    in_exponent = false
    for c in str
        k = findfirst(==(c), _SUPERSCRIPTS)
        if k === nothing
            in_exponent = false
            print(buf, c)
        else
            in_exponent || print(buf, '^')
            in_exponent = true
            print(buf, k - 1)
        end
    end
    str = String(take!(buf))
    rank = 0
    tors = BigInt[]
    for part in split(str, '+')
        if (m = match(r"^Z(?:\^(\d+))?$", part)) !== nothing
            rank += m[1] === nothing ? 1 : parse(Int, m[1])
        elseif (m = match(r"^Z/(\d+)Z?$", part)) !== nothing
            push!(tors, parse(BigInt, m[1]))
        elseif (m = match(r"^\(Z/(\d+)Z?\)\^(\d+)$", part)) !== nothing
            append!(tors, fill(parse(BigInt, m[1]), parse(Int, m[2])))
        elseif part != "0"
            throw(ArgumentError("cannot read \"$part\" in \"$s\": expected ℤ, ℤᵏ, ℤ/nℤ, (ℤ/nℤ)ᵏ or 0"))
        end
    end
    return FinitelyGeneratedAbelianGroup(rank, tors)
end
