# Tensor products of complexes over a monoid algebra by modules of coefficients, on one side
# or on both, and the homology with coefficients (J1.e, task e2).

"""
    GarsideHomology.tensor_product(N, C)
    GarsideHomology.tensor_product(C, M)
    GarsideHomology.tensor_product(N, C, M)
    N ⊗ C,  C ⊗ M                  (with `using GarsideHomology: ⊗`)

The tensor products over the monoid algebra of a complex `C` of free modules over a
[`MonoidAlgebra`](@ref) by a module of coefficients `N` on the left (as a right module) or
`M` on the right (as a left module), [`AbstractMonoidModule`](@ref)s over the ring of the
coefficients of `C` (ℤ for now).

- `C` of bimodules (an lcm resolution): `N ⊗ C` is the complex of free right modules whose
  terms `c·σ y τ` become `χ_N(σ)·c · y τ` (labels τ, side `:right`); `C ⊗ M` that of free
  left modules, `χ_M(τ)·c · σ y` (labels σ, side `:left`); in each column, the terms of the
  same row and label are summed. `tensor_product(N, C, M)`, or `(N ⊗ C) ⊗ M`, is the
  complex over the ring, `χ_N(σ)·χ_M(τ)·c · y`.
- `C` of right modules (`N ⊗ R`): `C ⊗ M` is the complex over the ring; of left modules:
  `N ⊗ C`. Tensoring a side that is already tensored is an `ArgumentError`.

The results share the bases of `C`, which keep their cells. `homology(R, M)` is the homology
of `TrivialModule ⊗ R ⊗ M`, that is `H₍*₎(B; M)` for the monoid `B` (the homology of its group
of fractions).
`tensor_product` is exported by Hecke and Oscar, `⊗` by Hecke: neither is exported by GarsideHomology.
"""
function tensor_product end

"""
    ⊗(N::AbstractMonoidModule, C), ⊗(C, M::AbstractMonoidModule)

The tensor products `GarsideHomology.tensor_product(N, C)` and `GarsideHomology.tensor_product(C, M)`. Public, not
exported: write `using GarsideHomology: ⊗`, which wins over the export of `⊗` by Hecke.
"""
⊗(N::AbstractMonoidModule, C::AbstractChainComplex) = tensor_product(N, C)
⊗(C::AbstractChainComplex, M::AbstractMonoidModule) = tensor_product(C, M)

function tensor_product(N::AbstractMonoidModule, C::AbstractChainComplex{<:MonoidAlgebra})
    A = base_ring(C)
    _check_coefficients(A, N)
    A.side === :right &&
        throw(ArgumentError("tensor_product: the left side of the complex is already tensored"))
    χ = _character(N, A.monoid)
    A.side === :left && return specialise(C, χ)
    return _one_sided(C, :right, l -> (χ(l[1]), l[2]))
end

function tensor_product(C::AbstractChainComplex{<:MonoidAlgebra}, M::AbstractMonoidModule)
    A = base_ring(C)
    _check_coefficients(A, M)
    A.side === :left &&
        throw(ArgumentError("tensor_product: the right side of the complex is already tensored"))
    χ = _character(M, A.monoid)
    A.side === :right && return specialise(C, χ)
    return _one_sided(C, :left, l -> (χ(l[2]), l[1]))
end

function tensor_product(N::AbstractMonoidModule, C::AbstractChainComplex{<:MonoidAlgebra},
                        M::AbstractMonoidModule)
    A = base_ring(C)
    A.side === :both || return tensor_product(tensor_product(N, C), M)
    _check_coefficients(A, N)
    _check_coefficients(A, M)
    χN, χM = _character(N, A.monoid), _character(M, A.monoid)
    return specialise(C, l -> χN(l[1]) * χM(l[2]))
end

tensor_product(::AbstractMonoidModule, ::AbstractChainComplex) =
    throw(ArgumentError("tensor_product: the complex is not over a monoid algebra"))
tensor_product(::AbstractChainComplex, ::AbstractMonoidModule) =
    throw(ArgumentError("tensor_product: the complex is not over a monoid algebra"))
tensor_product(::AbstractMonoidModule, ::AbstractChainComplex, ::AbstractMonoidModule) =
    throw(ArgumentError("tensor_product: the complex is not over a monoid algebra"))

# The module must be over the ring of the coefficients of the monoid algebra, the integers
# for now (other rings come with J2).
function _check_coefficients(A::MonoidAlgebra, M::AbstractMonoidModule)
    base_ring(M) isa AbstractAlgebra.Integers ||
        throw(ArgumentError("tensor_product: modules over the integers only for now, not over $(base_ring(M))"))
    base_ring(M) == base_ring(A) ||
        throw(ArgumentError("tensor_product: the module is over $(base_ring(M)), the complex over $(base_ring(A))"))
    rank(M) == 1 || throw(ArgumentError("tensor_product: modules of rank one only for now"))
    return nothing
end

# The complex of one-sided modules obtained from a complex of bimodules by f, which sends
# a label (σ, τ) to (weight, kept label): side :right keeps τ, side :left keeps σ.
function _one_sided(C::AbstractChainComplex{<:MonoidAlgebra}, side::Symbol, f)
    A = base_ring(C)
    B = MonoidAlgebra(base_ring(A), A.monoid, side)
    degrees = range(C)
    modules = [BasedFreeModule(B, cells(C[n])) for n in degrees]
    maps = [_one_sided_matrix(differential(C, n), f) for n in first(degrees):(last(degrees) - 1)]
    return ChainComplex(modules, maps; lo = first(degrees))
end

function _one_sided_matrix(d::LabelledSparseMatrix{NTuple{2, Int32}, Tv}, f) where {Tv}
    rows, labs, vals = rowvals(d), labels(d), nonzeros(d)
    I, J, L, V = Int[], Int[], Int32[], Tv[]
    for j in 1:size(d, 2), k in nzrange(d, j)
        w, l = f(labs[k])
        iszero(w) && continue
        push!(I, rows[k]); push!(J, j); push!(L, l)
        push!(V, Complexes._checked_mul(convert(Tv, w), vals[k]))
    end
    return labelled_sparse(I, J, L, V, size(d)...)
end

"""
    homology(R::AbstractChainComplex, M::AbstractMonoidModule; kwargs...)
    homology(R::AbstractChainComplex, M::AbstractMonoidModule, n; kwargs...)

The homology of the monoid, that of its group of fractions, with coefficients in the module
`M`, from a resolution `R` of bimodules (an lcm resolution): the homology over the ring of
`M` of `TrivialModule ⊗ R ⊗ M` ([`tensor_product`](@ref)), in all degrees or in degree `n`.
For a complex of right modules (`TrivialModule ⊗ R`), the homology of `C ⊗ M`. The keyword
arguments go to the solver, as for `homology(C)`.
"""
Complexes.homology(C::AbstractChainComplex{<:MonoidAlgebra}, M::AbstractMonoidModule; kwargs...) =
    homology(_with_coefficients(C, M); kwargs...)

Complexes.homology(C::AbstractChainComplex{<:MonoidAlgebra}, M::AbstractMonoidModule, n::Integer; kwargs...) =
    homology(_with_coefficients(C, M), n; kwargs...)

"""
    homology(G::AbstractGerm, M::AbstractMonoidModule = TrivialModule(); ordering = LCMOrdering(), kwargs...)
    homology(G::AbstractGerm, M::AbstractMonoidModule, n; ordering = LCMOrdering(), kwargs...)
    homology(G::AbstractGerm, n; ordering = LCMOrdering(), kwargs...)

The homology of the monoid of the germ `G`, that of its group of fractions, with
coefficients in the module `M` (trivial by default), in all degrees or in degree `n`: that
of its lcm resolution for `ordering`, `homology(lcm_resolution(G, ordering), M)`. The
resolution is built at each call: to compute several homologies of the same monoid, build
it once with [`lcm_resolution`](@ref). The other keyword arguments go to the solver, as
for `homology(C)`.

With Chevie loaded, the same forms for a `BraidMonoid` or a `DualBraidMonoid` of Chevie,
through its germ ([`germ`](@ref)): `homology(BraidMonoid(coxgroup(:A, 2)))` is the
homology of the braid group on three strands.
"""
Complexes.homology(G::AbstractGerm, M::AbstractMonoidModule = TrivialModule();
                   ordering::LCMOrdering = LCMOrdering(), kwargs...) =
    homology(lcm_resolution(G, ordering), M; kwargs...)

Complexes.homology(G::AbstractGerm, M::AbstractMonoidModule, n::Integer;
                   ordering::LCMOrdering = LCMOrdering(), kwargs...) =
    homology(lcm_resolution(G, ordering), M, n; kwargs...)

Complexes.homology(G::AbstractGerm, n::Integer; ordering::LCMOrdering = LCMOrdering(), kwargs...) =
    homology(lcm_resolution(G, ordering), TrivialModule(), n; kwargs...)

function _with_coefficients(C::AbstractChainComplex{<:MonoidAlgebra}, M::AbstractMonoidModule)
    A = base_ring(C)
    A.side === :left &&
        throw(ArgumentError("homology: a complex of left modules; tensor its left side, e.g. homology(TrivialModule(ZZ) ⊗ C)"))
    A.side === :right && return tensor_product(C, M)
    return tensor_product(TrivialModule(base_ring(M)), C, M)
end
