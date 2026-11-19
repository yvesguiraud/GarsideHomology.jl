# Tensor products of complexes over a monoid algebra by modules of coefficients, on one side
# or on both, and the homology with coefficients.

# A ring stands for its trivial module: `homology(R, ℤ)`, `ℤ ⊗ R ⊗ ℤ`.
const _Coefficients = Union{AbstractMonoidModule, AbstractAlgebra.Ring}
_module(M::AbstractMonoidModule) = M
_module(A::AbstractAlgebra.Ring) = TrivialModule(A)

"""
    GarsideHomology.tensor_product(N, C::AbstractChainComplex)
    GarsideHomology.tensor_product(C::AbstractChainComplex, M)
    GarsideHomology.tensor_product(N, C::AbstractChainComplex, M)
    N ⊗ C,  C ⊗ M                  (with `using GarsideHomology: ⊗`, or with Chevie loaded)

The tensor products over the monoid algebra of a complex `C` of free modules over a [`MonoidAlgebra`](@ref) by a module of coefficients `N` on the left (as a right module) or `M` on the right (as a left module). The arguments `N` and `M` are [`AbstractMonoidModule`](@ref)s over the ring of the coefficients of `C` (ℤ for now), or rings (`AbstractAlgebra.Ring`), a ring standing for its trivial module: `ℤ ⊗ C ⊗ ℤ` is `TrivialModule() ⊗ C ⊗ TrivialModule()`, that is, [`extend_scalars`](@ref)`(C)`.

The result depends on the sides of `C`:
- if `C` is a complex of bimodules, such as an lcm resolution, then `N ⊗ C` is the complex of free right modules whose terms `c·σ[y]τ` become `χ_N(σ)·c·[y]τ` (labels τ, side `:right`), and `C ⊗ M` is the complex of free left modules whose terms become `χ_M(τ)·c·σ[y]` (labels σ, side `:left`); in each column, the terms of the same row and label are summed; `tensor_product(N, C, M)`, or `(N ⊗ C) ⊗ M`, is the complex over the ring whose terms are `χ_N(σ)·χ_M(τ)·c·[y]`;
- if `C` is a complex of right modules, such as `N ⊗ R`, then `C ⊗ M` is the complex over the ring; if `C` is a complex of left modules, then `N ⊗ C` is the complex over the ring;
- tensoring a side that is already tensored throws an `ArgumentError`.

The results share the bases of `C`, which keep their cells. `homology(R, M)` is the homology of `TrivialModule() ⊗ R ⊗ M`, that is, the homology of the monoid with coefficients in `M` (see [`homology`](@ref)).

The name `tensor_product` is exported by Hecke and Oscar, and `⊗` by Hecke, Oscar and Chevie (that of FinitePosets), so that neither is exported by GarsideHomology. With Chevie loaded, the `⊗` of Chevie has the same methods for our complexes (package extension), so that `ℤ ⊗ R ⊗ ℤ` works without importing ours.
"""
function tensor_product end

"""
    ⊗(N, C::AbstractChainComplex)
    ⊗(C::AbstractChainComplex, M)

The tensor products `GarsideHomology.tensor_product(N, C)` and `GarsideHomology.tensor_product(C, M)`, where `N` and `M` are modules of coefficients or rings (for their trivial modules). The name is public but not exported: write `using GarsideHomology: ⊗`, which wins over the exports of `⊗` by Hecke, Oscar and Chevie. With Chevie loaded, the `⊗` of Chevie (that of FinitePosets) has the same methods for our complexes, and `ℤ ⊗ R ⊗ ℤ` works without the import.
"""
⊗(N::_Coefficients, C::AbstractChainComplex) = tensor_product(N, C)
⊗(C::AbstractChainComplex, M::_Coefficients) = tensor_product(C, M)

# The arguments are not reassigned (a module for a ring): the local names left and right
# keep a single type each.
# Signature line, counted only by some compiled instances in the coverage of Julia.
function tensor_product(N::_Coefficients, C::AbstractChainComplex{<:MonoidAlgebra})  # COV_EXCL_LINE
    left = _module(N)
    A = base_ring(C)
    _check_coefficients(A, left)
    A.side === :right &&
        throw(ArgumentError("tensor_product: the left side of the complex is already tensored"))
    χ = _character(left, A.monoid)
    A.side === :left && return extend_scalars(C, χ)
    return _one_sided(C, :right, l -> (χ(l[1]), l[2]))
end

# Signature line, counted only by some compiled instances in the coverage of Julia.
function tensor_product(C::AbstractChainComplex{<:MonoidAlgebra}, M::_Coefficients)  # COV_EXCL_LINE
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
# for now (other rings will come later).
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
    homology(R::AbstractChainComplex{<:MonoidAlgebra}, M; kwargs...)
    homology(R::AbstractChainComplex{<:MonoidAlgebra}, M, n::Integer; kwargs...)

The homology of the monoid with coefficients in the module `M`, computed from a resolution `R` of bimodules, such as an lcm resolution, in all degrees or in degree `n`. It is the homology over the ring of `M` of `TrivialModule() ⊗ R ⊗ M` ([`tensor_product`](@ref)). Since the monoid is a Garside monoid, this is also the homology of its enveloping group with coefficients in `M` (see the page Mathematical background of the User guide).

The argument `M` is a module of coefficients ([`AbstractMonoidModule`](@ref)), or a ring for its trivial module: `homology(R, ℤ)` is `homology(R, TrivialModule())`. For a complex of right modules, such as `TrivialModule() ⊗ R`, the result is the homology of `C ⊗ M`. The keyword arguments go to the solver, as for `homology(C)`.

Without coefficients, `homology(R)` is an error: the homology of the resolution itself is not that of the monoid.
"""
Complexes.homology(C::AbstractChainComplex{<:MonoidAlgebra}, M::_Coefficients; kwargs...) =
    homology(_with_coefficients(C, _module(M)); kwargs...)

Complexes.homology(C::AbstractChainComplex{<:MonoidAlgebra}, M::_Coefficients, n::Integer; kwargs...) =
    homology(_with_coefficients(C, _module(M)), n; kwargs...)

"""
    homology(G::AbstractGerm, M = TrivialModule(); ordering::LCMOrdering = LCMOrdering(), kwargs...)
    homology(G::AbstractGerm, M, n::Integer; ordering::LCMOrdering = LCMOrdering(), kwargs...)
    homology(G::AbstractGerm, n::Integer; ordering::LCMOrdering = LCMOrdering(), kwargs...)

The homology of the monoid of the germ `G` with coefficients in the module `M`, in all degrees or in degree `n`. The module is trivial by default, and a ring stands for its trivial module, as in `homology(G, ℤ)`. Since the monoid is a Garside monoid, this is also the homology of its enveloping group.

The homology is computed from the lcm resolution for `ordering`, as `homology(lcm_resolution(G, ordering), M)`. The resolution is built at each call: to compute several homologies of the same monoid, build it once with [`lcm_resolution`](@ref). The other keyword arguments go to the solver, as for `homology(C)`.

With Chevie loaded, the same forms apply to a `BraidMonoid` or a `DualBraidMonoid` of Chevie, through its germ ([`germ`](@ref)): for instance, `homology(BraidMonoid(coxgroup(:A, 2)))` is the homology of the braid group on three strands.
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
