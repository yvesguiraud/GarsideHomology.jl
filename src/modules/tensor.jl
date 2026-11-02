# Tensor products of complexes over a monoid algebra by modules of coefficients, on one side
# or on both, and the homology with coefficients (J1.e, task e2).

# A ring stands for its trivial module (J1.g): `homology(R, ℤ)`, `ℤ ⊗ R ⊗ ℤ`.
const _Coefficients = Union{AbstractMonoidModule, AbstractAlgebra.Ring}
_module(M::AbstractMonoidModule) = M
_module(A::AbstractAlgebra.Ring) = TrivialModule(A)

"""
    GarsideHomology.tensor_product(N, C)
    GarsideHomology.tensor_product(C, M)
    GarsideHomology.tensor_product(N, C, M)
    N ⊗ C,  C ⊗ M                  (with `using GarsideHomology: ⊗`, or with Chevie loaded)

The tensor products over the monoid algebra of a complex `C` of free modules over a
[`MonoidAlgebra`](@ref) by a module of coefficients `N` on the left (as a right module) or
`M` on the right (as a left module), [`AbstractMonoidModule`](@ref)s over the ring of the
coefficients of `C` (ℤ for now). A ring stands for its trivial module: `ℤ ⊗ C ⊗ ℤ` is
`TrivialModule() ⊗ C ⊗ TrivialModule()`, that is [`extend_scalars`](@ref)`(C)`.

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
`tensor_product` is exported by Hecke and Oscar, `⊗` by Hecke, Oscar and Chevie (that of
FinitePosets): neither is exported by GarsideHomology. With Chevie loaded, the `⊗` of
Chevie has the same methods for our complexes (package extension), so that `ℤ ⊗ R ⊗ ℤ`
works without importing ours.
"""
function tensor_product end

"""
    ⊗(N, C), ⊗(C, M)

The tensor products `GarsideHomology.tensor_product(N, C)` and
`GarsideHomology.tensor_product(C, M)`, `N` and `M` modules of coefficients or rings (for
their trivial modules). Public, not exported: write `using GarsideHomology: ⊗`, which wins
over the exports of `⊗` by Hecke, Oscar and Chevie. With Chevie loaded, the `⊗` of Chevie
(that of FinitePosets) has the same methods for our complexes, and `ℤ ⊗ R ⊗ ℤ` works
without the import.
"""
⊗(N::_Coefficients, C::AbstractChainComplex) = tensor_product(N, C)
⊗(C::AbstractChainComplex, M::_Coefficients) = tensor_product(C, M)

# The arguments are not reassigned (a module for a ring): the local names left and right
# keep a single type each.
function tensor_product(N::_Coefficients, C::AbstractChainComplex{<:MonoidAlgebra})
    left = _module(N)
    A = base_ring(C)
    _check_coefficients(A, left)
    A.side === :right &&
        throw(ArgumentError("tensor_product: the left side of the complex is already tensored"))
    χ = _character(left, A.monoid)
    A.side === :left && return extend_scalars(C, χ)
    return _one_sided(C, :right, l -> (χ(l[1]), l[2]))
end

function tensor_product(C::AbstractChainComplex{<:MonoidAlgebra}, M::_Coefficients)
    right = _module(M)
    A = base_ring(C)
    _check_coefficients(A, right)
    A.side === :left &&
        throw(ArgumentError("tensor_product: the right side of the complex is already tensored"))
    χ = _character(right, A.monoid)
    A.side === :right && return extend_scalars(C, χ)
    return _one_sided(C, :left, l -> (χ(l[2]), l[1]))
end

function tensor_product(N::_Coefficients, C::AbstractChainComplex{<:MonoidAlgebra},
                        M::_Coefficients)
    left, right = _module(N), _module(M)
    A = base_ring(C)
    A.side === :both || return tensor_product(tensor_product(left, C), right)
    _check_coefficients(A, left)
    _check_coefficients(A, right)
    χN, χM = _character(left, A.monoid), _character(right, A.monoid)
    return extend_scalars(C, l -> χN(l[1]) * χM(l[2]))
end

tensor_product(::_Coefficients, ::AbstractChainComplex) =
    throw(ArgumentError("tensor_product: the complex is not over a monoid algebra"))
tensor_product(::AbstractChainComplex, ::_Coefficients) =
    throw(ArgumentError("tensor_product: the complex is not over a monoid algebra"))
tensor_product(::_Coefficients, ::AbstractChainComplex, ::_Coefficients) =
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
    modules = [BasedFreeModule(B, basis(C[n])) for n in degrees]
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
    return LabelledSparseMatrix(I, J, L, V, size(d)...)
end

"""
    homology(R::AbstractChainComplex, M; kwargs...)
    homology(R::AbstractChainComplex, M, n; kwargs...)

The homology of the monoid, that of its group of fractions, with coefficients in the module
`M`, from a resolution `R` of bimodules (an lcm resolution): the homology over the ring of
`M` of `TrivialModule ⊗ R ⊗ M` ([`tensor_product`](@ref)), in all degrees or in degree `n`.
For a complex of right modules (`TrivialModule ⊗ R`), the homology of `C ⊗ M`. `M` is a
module of coefficients, or a ring for its trivial module: `homology(R, ℤ)` is
`homology(R, TrivialModule())`. The keyword arguments go to the solver, as for
`homology(C)`. Without coefficients, `homology(R)` is an error: the homology of the
resolution itself is not that of the monoid.
"""
Complexes.homology(C::AbstractChainComplex{<:MonoidAlgebra}, M::_Coefficients; kwargs...) =
    homology(_with_coefficients(C, _module(M)); kwargs...)

Complexes.homology(C::AbstractChainComplex{<:MonoidAlgebra}, M::_Coefficients, n::Integer; kwargs...) =
    homology(_with_coefficients(C, _module(M)), n; kwargs...)

"""
    homology(G::AbstractGerm, M = TrivialModule(); ordering = LCMOrdering(), kwargs...)
    homology(G::AbstractGerm, M, n; ordering = LCMOrdering(), kwargs...)
    homology(G::AbstractGerm, n; ordering = LCMOrdering(), kwargs...)

The homology of the monoid of the germ `G`, that of its group of fractions, with
coefficients in the module `M` (trivial by default; a ring stands for its trivial module,
as in `homology(G, ℤ)`), in all degrees or in degree `n`: that
of its lcm resolution for `ordering`, `homology(lcm_resolution(G, ordering), M)`. The
resolution is built at each call: to compute several homologies of the same monoid, build
it once with [`lcm_resolution`](@ref). The other keyword arguments go to the solver, as
for `homology(C)`.

With Chevie loaded, the same forms for a `BraidMonoid` or a `DualBraidMonoid` of Chevie,
through its germ ([`germ`](@ref)): `homology(BraidMonoid(coxgroup(:A, 2)))` is the
homology of the braid group on three strands.
"""
Complexes.homology(G::AbstractGerm, M::_Coefficients = TrivialModule();
                   ordering::LCMOrdering = LCMOrdering(), kwargs...) =
    homology(lcm_resolution(G, ordering), M; kwargs...)

Complexes.homology(G::AbstractGerm, M::_Coefficients, n::Integer;
                   ordering::LCMOrdering = LCMOrdering(), kwargs...) =
    homology(lcm_resolution(G, ordering), M, n; kwargs...)

Complexes.homology(G::AbstractGerm, n::Integer; ordering::LCMOrdering = LCMOrdering(), kwargs...) =
    homology(lcm_resolution(G, ordering), TrivialModule(), n; kwargs...)

function _with_coefficients(C::AbstractChainComplex{<:MonoidAlgebra}, M::AbstractMonoidModule)
    A = base_ring(C)
    A.side === :left &&
        throw(ArgumentError("homology: a complex of left modules; tensor its left side, e.g. homology(ℤ ⊗ C)"))
    A.side === :right && return tensor_product(C, M)
    return tensor_product(TrivialModule(base_ring(M)), C, M)
end
