# Certificates of the lcm resolutions (J1.e, task e1): homogeneity of the terms, ∂² = 0 in
# the bimodules and Euler characteristic by weight (correctness); no term 1·y·1, cells of
# each weight in a single degree, lengths of the labels (properties).

"""
    ResolutionCertificate

The report of [`certify`](@ref) on a resolution of bimodules over a germ. Each check gives
a number of failures, `0` when it passes, or `nothing` when it was not done.

Correctness, which must hold for every ordering, minimal or not:
- `homogeneity`: the terms `σ y τ` of a boundary `∂x` with `σ·π(y)·τ ≠ π(x)`, `π` being the
  weight of a cell;
- `square_zero`: the columns of `d_n ∘ d_{n+1}` that are not zero, computed in the
  bimodules with the products of the germ (`nothing` with `square_zero = false`);
- `euler`: the simples `ω` where the Euler characteristic of the cells of weight `ω`,
  `Σₙ (−1)ⁿ cₙ(ω)`, differs from `μ(1, ω)` (`nothing` without the Möbius function).

Properties, never failures of the resolution:
- `identity_terms[n + 1]`: the terms `1·y·1` of `d_n`, which a minimal resolution has none;
- `split_weights`: the weights whose cells lie in several degrees (`nothing` without μ);
  with the Euler characteristic, none means `cₙ(ω) = |μ(1, ω)|` in a single degree;
- `mobius_sum`: `Σ_ω |μ(1, ω)|`, the number of cells of such a resolution;
- `max_label_length[n + 1]`: the largest `ℓ(σ) + ℓ(τ)` over the terms of `d_n` (0 if there
  are none).

Also: `cells[n + 1]` and `terms[n + 1]`, the numbers of cells of degree `n` and of terms of
`d_n`; `mobius_method`, how μ was computed (`:crosscut`, `:layers` or `:none`);
`failures`, the first failure of each failed check, in words (degree, cell named by its
atom sequence as in the snapshots, term or weight); `times`, the seconds of each step.

Read it with [`is_correct`](@ref), [`is_minimal`](@ref) and [`is_linear`](@ref).
"""
struct ResolutionCertificate
    cells::Vector{Int}
    terms::Vector{Int}
    homogeneity::Int
    square_zero::Union{Int, Nothing}
    mobius_method::Symbol
    euler::Union{Int, Nothing}
    mobius_sum::Union{Int, Nothing}
    split_weights::Union{Int, Nothing}
    identity_terms::Vector{Int}
    max_label_length::Vector{Int}
    failures::Vector{Pair{Symbol, String}}
    times::Vector{Pair{Symbol, Float64}}
end

"""
    is_correct(c::ResolutionCertificate) -> Bool

Whether the three checks of correctness were done and passed: homogeneity, ∂² = 0 and the
Euler characteristic by weight. A check that was not done gives `false`: the resolution is
then not certified.
"""
is_correct(c::ResolutionCertificate) = c.homogeneity == 0 && c.square_zero === 0 && c.euler === 0

"""
    is_minimal(c::ResolutionCertificate) -> Bool

Whether the differential has no term `1·y·1`: the criterion of minimality of a resolution
graded by the weight.
"""
is_minimal(c::ResolutionCertificate) = all(iszero, c.identity_terms)

"""
    is_linear(c::ResolutionCertificate) -> Bool

Whether every term `σ y τ` of the differential has `ℓ(σ) + ℓ(τ) = 1`: no term `1·y·1` and
at most 1 in each degree.
"""
is_linear(c::ResolutionCertificate) = is_minimal(c) && all(<=(1), c.max_label_length)

# Above this number of atoms, mobius = :auto does not use the crosscut (2^N sets of atoms);
# above the second one, mobius = :crosscut is refused.
const _CROSSCUT_AUTO = 16
const _CROSSCUT_MAX = 24

"""
    certify(R::LCMResolution; square_zero = true, mobius = :auto,
            nworkers = Threads.nthreads()) -> ResolutionCertificate

Check the lcm resolution `R` ([`lcm_resolution`](@ref)) and return the report, a
[`ResolutionCertificate`](@ref); a failed check is a field of the report, never an
exception.

- Homogeneity, terms `1·y·1` and lengths of the labels: one pass over the terms.
- ∂² = 0 in the bimodules: each column of `d_{n+1}` composed with the columns of `d_n`, the
  products of simples read in the germ, then put in normal form. With `nworkers > 1` and a
  germ safe for concurrent reads (a `TableGerm`), the columns are split between `nworkers`
  tasks; a lazy germ stays on one thread. `square_zero = false` skips this check, the
  costliest.
- The Möbius function μ(1, ω), for the Euler characteristic and the degrees by weight:
  `mobius = :crosscut` by the crosscut theorem of Rota (2^N lcms, any germ, at most 24
  atoms), `:layers` by a pass over the intervals (a germ that lists its simples), `:none`
  to skip; `:auto` takes the crosscut up to 16 atoms (the braid monoids of Coxeter groups,
  even lazy), the layers above (the dual braid monoids).

The arithmetic of the coefficients is checked: an overflow throws an `OverflowError`
rather than a wrong report.
"""
function certify(R::LCMResolution; square_zero::Bool = true, mobius::Symbol = :auto,
                 nworkers::Integer = Threads.nthreads())
    nworkers >= 1 || throw(ArgumentError("certify: nworkers = $nworkers, at least 1 expected"))
    T = cell_tree(R)
    G = T.germ
    method = _mobius_method(G, mobius)
    top = last(range(R))
    failures = Pair{Symbol, String}[]
    times = Pair{Symbol, Float64}[]
    cells = [ncells(T, n) for n in 0:top]
    terms = [nnz(differential(R, n)) for n in 0:(top - 1)]

    t = time_ns()
    homogeneity, identity_terms, max_length, h_failure = _homogeneity(R, T)
    h_failure === nothing || push!(failures, :homogeneity => h_failure)
    push!(times, :homogeneity => _seconds(t))

    square = nothing
    if square_zero
        t = time_ns()
        square, s_failure = _square_zero(R, G; nworkers)
        s_failure === nothing || push!(failures, :square_zero => _square_zero_message(T, s_failure))
        push!(times, :square_zero => _seconds(t))
    end

    euler = mobius_sum = split_weights = nothing
    if method !== :none
        t = time_ns()
        μ = method === :crosscut ? _mobius_crosscut(G) : _mobius_values(G)
        push!(times, :mobius => _seconds(t))
        t = time_ns()
        euler, mobius_sum, split_weights, e_failure = _euler(T, μ)
        e_failure === nothing || push!(failures, :euler => e_failure)
        push!(times, :euler => _seconds(t))
    end
    return ResolutionCertificate(cells, terms, homogeneity, square, method, euler, mobius_sum, split_weights,
                                 identity_terms, max_length, failures, times)
end

certify(::AbstractChainComplex; kwargs...) =
    throw(ArgumentError("certify: the complex is not an lcm resolution (see lcm_resolution)"))

function _mobius_method(G::AbstractGerm, mobius::Symbol)
    N = natoms(G)
    if mobius === :auto
        return N <= _CROSSCUT_AUTO ? :crosscut : applicable(simples, G) ? :layers : :none
    elseif mobius === :crosscut
        N <= _CROSSCUT_MAX ||
            throw(ArgumentError("certify: mobius = :crosscut needs the 2^$N sets of atoms; at most $_CROSSCUT_MAX atoms"))
    elseif mobius === :layers
        applicable(simples, G) ||
            throw(ArgumentError("certify: mobius = :layers needs a germ that lists its simples (a tabulated germ)"))
    elseif mobius !== :none
        throw(ArgumentError("certify: mobius = :$mobius, expected :auto, :crosscut, :layers or :none"))
    end
    return mobius
end

_seconds(t::UInt64) = (time_ns() - t) / 1e9

# The product of two simples, the identity 0 without a call to the germ; nothing if it is
# not simple.
@inline function _mul_or_nothing(G::AbstractGerm, s::Integer, t::Integer)
    iszero(s) && return Int32(t)
    iszero(t) && return Int32(s)
    return product(G, s, t)
end

# Names in the messages: a cell by its atom sequence ("-" for the cell of degree 0), as in
# the snapshots; a term as in their D lines.
_cell_string(T::CellTree, x::Integer) = cell_degree(T, x) == 0 ? "-" : join(atom_sequence(T, x), '.')
_term_string(G::AbstractGerm, T::CellTree, c::Integer, σ::Integer, y::Integer, τ::Integer) =
    string(c > 0 ? "+$c" : "$c", " ", _simple_string(G, σ), " [", _cell_string(T, y), "] ", _simple_string(G, τ))

# ------------------------------------------------------------------------------------
# Homogeneity, terms 1·y·1, lengths of the labels

# One pass over the terms of the differential: the number of terms σ y τ of ∂x with
# σ·π(y)·τ ≠ π(x), the terms 1·y·1 and the largest ℓ(σ) + ℓ(τ) of each d_n, the first
# inhomogeneous term in words (or nothing).
function _homogeneity(R::LCMResolution, T::CellTree)
    G = T.germ
    top = last(range(R))
    bad = 0
    first_bad = nothing
    identity = zeros(Int, top)
    max_length = zeros(Int, top)
    for n in 0:(top - 1)
        d = differential(R, n)
        rows, labs, vals = rowvals(d), labels(d), nonzeros(d)
        row0 = first(cell_range(T, n)) - 1
        col0 = first(cell_range(T, n + 1)) - 1
        for j in 1:size(d, 2)
            x = col0 + j
            πx = cell_weight(T, x)
            for k in nzrange(d, j)
                y = row0 + Int(rows[k])
                σ, τ = labs[k]
                max_length[n + 1] = max(max_length[n + 1], length(G, σ) + length(G, τ))
                iszero(σ) && iszero(τ) && (identity[n + 1] += 1)
                p = _mul_or_nothing(G, σ, cell_weight(T, y))
                p === nothing || (p = _mul_or_nothing(G, p, τ))
                if p != πx
                    bad += 1
                    first_bad === nothing && (first_bad = string("d_$n, cell ", _cell_string(T, x), ": term ",
                                                                 _term_string(G, T, vals[k], σ, y, τ),
                                                                 " of another weight"))
                end
            end
        end
    end
    return bad, identity, max_length, first_bad
end

# ------------------------------------------------------------------------------------
# ∂² = 0 in the bimodules

# The columns of d_n ∘ d_{n+1} that are not zero, over all n, for a complex of bimodules
# over the germ G: their number and the first one in the order (n, column), as
# (n, column, info), info being its first remaining term (a BimoduleTerm whose cell is a
# row of d_n) or the pair of simples whose product is not simple; nothing if there is none.
# The columns are composed independently, split between nworkers tasks when the germ is
# safe for concurrent reads; the result does not depend on the split.
function _square_zero(C::AbstractChainComplex, G::AbstractGerm; nworkers::Integer = 1)
    nw = _thread_safe(G) ? Int(nworkers) : 1
    degrees = range(C)
    total = 0
    first_failure = nothing
    for n in first(degrees):(last(degrees) - 2)
        count, j, info = _square_zero_pair(differential(C, n), differential(C, n + 1), G, nw)
        total += count
        first_failure === nothing && j > 0 && (first_failure = (n = n, column = j, info = info))
    end
    return total, first_failure
end

function _square_zero_pair(A::LabelledSparseMatrix, B::LabelledSparseMatrix, G::AbstractGerm, nw::Int)
    ncol = size(B, 2)
    (nw == 1 || ncol <= 1) && return _square_zero_columns(A, B, G, 1:ncol, BimoduleTerm[])
    chunks = collect(Iterators.partition(1:ncol, cld(ncol, 8 * nw)))
    next = Threads.Atomic{Int}(1)
    tasks = [Threads.@spawn(_square_zero_worker(A, B, G, chunks, next)) for _ in 1:min(nw, length(chunks))]
    count, first_j, first_info = 0, 0, nothing
    for task in tasks
        c, j, info = fetch(task)
        count += c
        if j > 0 && (first_j == 0 || j < first_j)
            first_j, first_info = j, info
        end
    end
    return count, first_j, first_info
end

# A task: the chunks of columns taken in turn from the shared counter.
function _square_zero_worker(A, B, G, chunks, next)
    buffer = BimoduleTerm[]
    count, first_j, first_info = 0, 0, nothing
    while true
        i = Threads.atomic_add!(next, 1)
        i > length(chunks) && break
        c, j, info = _square_zero_columns(A, B, G, chunks[i], buffer)
        count += c
        if j > 0 && (first_j == 0 || j < first_j)
            first_j, first_info = j, info
        end
    end
    return count, first_j, first_info
end

# The columns cols of d_n ∘ d_{n+1}, in increasing order: the number of those that are not
# zero, the first one and what remains in it.
function _square_zero_columns(A, B, G, cols, buffer)
    count, first_j, first_info = 0, 0, nothing
    for j in cols
        info = _compose_column!(buffer, A, B, j, G)
        info === nothing && continue
        count += 1
        if first_j == 0
            first_j, first_info = j, info
        end
    end
    return count, first_j, first_info
end

# The column j of A ∘ B in normal form in buffer: c·c′·(σσ′) z (τ′τ) for each term c σ y τ
# of the column j of B and each term c′ σ′ z τ′ of the column y of A. Returns nothing if it
# is zero, its first term otherwise, or the pair of simples whose product is not simple.
function _compose_column!(buffer::Vector{BimoduleTerm}, A::LabelledSparseMatrix,
                          B::LabelledSparseMatrix, j::Integer, G::AbstractGerm)
    empty!(buffer)
    Ar, Al, Av = rowvals(A), labels(A), nonzeros(A)
    Br, Bl, Bv = rowvals(B), labels(B), nonzeros(B)
    for q in nzrange(B, j)
        y = Br[q]
        σ, τ = Bl[q]
        c = Int64(Bv[q])
        for p in nzrange(A, y)
            σ′, τ′ = Al[p]
            s = _mul_or_nothing(G, σ, σ′)
            s === nothing && return (Int32(σ), Int32(σ′))
            t = _mul_or_nothing(G, τ′, τ)
            t === nothing && return (Int32(τ′), Int32(τ))
            push!(buffer, BimoduleTerm(Ar[p], s, t, Base.checked_mul(c, Int64(Av[p]))))
        end
    end
    _normalise!(buffer)
    return isempty(buffer) ? nothing : buffer[1]
end

function _square_zero_message(T::CellTree, f)
    G = T.germ
    x = first(cell_range(T, f.n + 2)) + f.column - 1
    head = string("d_", f.n, " ∘ d_", f.n + 1, ", cell ", _cell_string(T, x), ": ")
    info = f.info
    info isa BimoduleTerm || return string(head, "the product of ", _simple_string(G, info[1]), " and ",
                                           _simple_string(G, info[2]), " is not simple")
    z = first(cell_range(T, f.n)) + Int(info.cell) - 1
    return string(head, "remaining term ", _term_string(G, T, info.coef, info.left, z, info.right))
end

"""
    is_square_zero(C::AbstractChainComplex; nworkers = Threads.nthreads())

For a complex of bimodules over a germ (an lcm resolution, for instance): whether
`d_n ∘ d_{n+1} = 0` for every `n`, computed in the bimodules with the products of the germ,
as by [`certify`](@ref), with `nworkers` tasks for a germ safe for concurrent reads.
"""
function Complexes.is_square_zero(C::AbstractChainComplex{A, _LCMMatrix};
                                  nworkers::Integer = Threads.nthreads()) where {A <: MonoidAlgebra{<:Any, <:AbstractGerm}}
    algebra = base_ring(C)
    algebra.side === :both || throw(ArgumentError("is_square_zero: labels in pairs, but the algebra is one-sided"))
    nworkers >= 1 || throw(ArgumentError("is_square_zero: nworkers = $nworkers, at least 1 expected"))
    return first(_square_zero(C, algebra.monoid; nworkers)) == 0
end

# ------------------------------------------------------------------------------------
# Euler characteristic and degrees by weight

_mu(μ::Vector{Int}, ω::Integer) = μ[Int(ω) + 1]
_mu(μ::Dict{Int32, Int}, ω::Integer) = get(μ, Int32(ω), 0)
_mu_support(μ::Vector{Int}) = (Int32(i - 1) for i in eachindex(μ) if !iszero(μ[i]))
_mu_support(μ::Dict{Int32, Int}) = keys(μ)
_mu_abs_sum(μ::Vector{Int}) = sum(abs, μ)
_mu_abs_sum(μ::Dict{Int32, Int}) = sum(abs, values(μ); init = 0)

# The weights where Σₙ (−1)ⁿ cₙ(ω) ≠ μ(1, ω), among those of the cells and those where μ is
# not zero; Σ |μ(1, ω)|; the weights whose cells lie in several degrees; the first failure
# in words (or nothing). μ is the vector of all the values or the dictionary of the nonzero
# ones.
function _euler(T::CellTree, μ)
    G = T.germ
    stats = Dict{Int32, NTuple{3, Int}}()    # weight => (Σ (−1)ⁿ cₙ, lowest degree, highest degree)
    for x in 1:ncells(T)
        ω = cell_weight(T, x)
        d = cell_degree(T, x)
        e, lo, hi = get(stats, ω, (0, typemax(Int), -1))
        stats[ω] = (iseven(d) ? e + 1 : e - 1, min(lo, d), max(hi, d))
    end
    split_weights = count(s -> s[2] != s[3], values(stats))
    weights = Set{Int32}(keys(stats))
    union!(weights, _mu_support(μ))
    bad = sort!([ω for ω in weights if first(get(stats, ω, (0, 0, 0))) != _mu(μ, ω)])
    message = nothing
    if !isempty(bad)
        ω = first(bad)
        message = string("weight ", _simple_string(G, ω), ": Σₙ (−1)ⁿ cₙ = ", first(get(stats, ω, (0, 0, 0))),
                         ", μ(1, ω) = ", _mu(μ, ω))
    end
    return length(bad), _mu_abs_sum(μ), split_weights, message
end

# ------------------------------------------------------------------------------------
# Display

_plural(k::Integer, word::AbstractString) = k == 1 ? "1 $word" : "$k $(word)s"

function _correctness_summary(c::ResolutionCertificate)
    is_correct(c) && return "correct"
    parts = String[]
    c.homogeneity > 0 && push!(parts, _plural(c.homogeneity, "inhomogeneous term"))
    if c.square_zero === nothing
        push!(parts, "∂² not checked")
    elseif c.square_zero > 0
        push!(parts, "∂² ≠ 0 on " * _plural(c.square_zero, "column"))
    end
    if c.euler === nothing
        push!(parts, "Euler not checked")
    elseif c.euler > 0
        push!(parts, "Euler fails on " * _plural(c.euler, "weight"))
    end
    all(endswith(" not checked"), parts) && return "not certified (" * join(parts, ", ") * ")"
    return "not correct (" * join(parts, ", ") * ")"
end

_minimality_summary(c::ResolutionCertificate) =
    is_minimal(c) ? "minimal" : "not minimal (" * _plural(sum(c.identity_terms), "term") * " 1·y·1)"

function _linearity_summary(c::ResolutionCertificate)
    is_linear(c) && return "linear"
    m = maximum(c.max_label_length; init = 0)
    return m <= 1 ? "not linear (terms 1·y·1)" : "not linear (max ℓ(σ) + ℓ(τ) = $m)"
end

Base.show(io::IO, c::ResolutionCertificate) =
    print(io, "ResolutionCertificate: ", _correctness_summary(c), ", ", _minimality_summary(c), ", ",
          _linearity_summary(c))

_time(c::ResolutionCertificate, steps...) =
    any(k -> any(p -> p.first === k, c.times), steps) ?
        sum((p.second for p in c.times if p.first in steps); init = 0.0) : nothing

_failure(c::ResolutionCertificate, key::Symbol) =
    (i = findfirst(p -> p.first === key, c.failures); i === nothing ? nothing : c.failures[i].second)

function _check_line(io::IO, c::ResolutionCertificate, key::Symbol, label::String, status::String, steps...)
    t = _time(c, key, steps...)
    print(io, "\n    ", rpad(label, 24), rpad(status, 24), t === nothing ? "" : string(round(t; digits = 2), " s"))
    message = _failure(c, key)
    message === nothing || print(io, "\n      first: ", message)
end

function Base.show(io::IO, ::MIME"text/plain", c::ResolutionCertificate)
    print(io, "ResolutionCertificate: ", sum(c.cells), " cells, ", sum(c.terms; init = 0), " terms")
    print(io, "\n  correct: ", is_correct(c) ? "yes" : "no: " * _correctness_summary(c))
    _check_line(io, c, :homogeneity, "homogeneity",
                c.homogeneity == 0 ? "ok" : _plural(c.homogeneity, "term"))
    _check_line(io, c, :square_zero, "∂² = 0 (bimodules)",
                c.square_zero === nothing ? "not checked" : c.square_zero == 0 ? "ok" : _plural(c.square_zero, "column"))
    _check_line(io, c, :euler, "Euler by weight",
                c.euler === nothing ? "not checked" :
                    (c.euler == 0 ? "ok" : _plural(c.euler, "weight")) * " (μ by $(c.mobius_method))", :mobius)
    print(io, "\n  minimal: ", is_minimal(c) ? "yes (no term 1·y·1)" :
              "no (" * _plural(sum(c.identity_terms), "term") * " 1·y·1, by degree " * join(c.identity_terms, ", ") * ")")
    print(io, "\n    ", rpad("one degree per weight", 24))
    if c.split_weights === nothing
        print(io, "not checked")
    else
        print(io, c.split_weights == 0 ? "yes" : "no, " * _plural(c.split_weights, "weight") * " in several degrees",
              " (Σ|μ| = ", c.mobius_sum, ", ", sum(c.cells), " cells)")
        is_minimal(c) && c.split_weights > 0 && c.euler == 0 &&
            print(io, "\n    note: minimal, but some intervals have homology in several degrees")
    end
    print(io, "\n  linear: ", is_linear(c) ? "yes" : "no", " (max ℓ(σ) + ℓ(τ) by degree: ",
          isempty(c.max_label_length) ? "none" : join(c.max_label_length, ", "), ")")
end
