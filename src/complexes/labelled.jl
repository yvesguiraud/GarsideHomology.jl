# Labelled sparse matrices: sparse matrices over a monoid algebra.

"""
    LabelledSparseMatrix{L, Tv}

Sparse `m × n` matrix over a monoid algebra ([`MonoidAlgebra`](@ref)). Writing the
matrix as `D = Σ_γ D_γ · γ`, with integer matrices `D_γ` and `γ` in the monoid, all the
`D_γ` are stored in a single CSC structure, each stored integer carrying its `γ`, its
label. Column `j` is the image of the `j`-th cell of the source; its stored entry `k` is
the monomial `nzval[k] · label[k]` in row `rowval[k]`, and the coefficient in row `i` is
the sum of the monomials of column `j` in row `i`. A label is the code of an element of
the monoid (`Int32`) or of a pair `(σ, τ)` for a bimodule (`NTuple{2, Int32}`), so that a
term of a bimodule differential reads `c · σ y τ`.

Invariants, checked by the constructor unless `check = false` (see
[`check_invariants`](@ref)):

1. `length(colptr) == n + 1`, `colptr[1] == 1`, `colptr` is nondecreasing and
   `colptr[n + 1] - 1` is the common length of `rowval`, `label` and `nzval`;
2. `1 ≤ rowval[k] ≤ m`;
3. in each column, the pairs `(rowval[k], label[k])` are strictly increasing: the
   storage is canonical, two equal matrices have equal arrays;
4. no stored coefficient is zero;
5. labels are nonnegative, label `0` (or `(0, 0)`) being the identity.

Constructors: `LabelledSparseMatrix(m, n, colptr, rowval, label, nzval; check = true)`
from the CSC arrays, `LabelledSparseMatrix{L, Tv}(m, n)` for the zero matrix, and
[`labelled_sparse`](@ref) from a list of monomials. Accessors follow SparseArrays:
`size`, `nrows`, `ncols`, `nnz`, `nzrange`, `rowvals`, `nonzeros`, plus [`labels`](@ref).
"""
struct LabelledSparseMatrix{L <: Union{Integer, NTuple{2, Integer}}, Tv}
    m::Int
    n::Int
    colptr::Vector{Int}
    rowval::Vector{Int32}
    label::Vector{L}
    nzval::Vector{Tv}
    function LabelledSparseMatrix{L, Tv}(m::Integer, n::Integer, colptr::Vector{Int},
                                         rowval::Vector{Int32}, label::Vector{L},
                                         nzval::Vector{Tv}; check::Bool = true) where {L, Tv}
        d = new{L, Tv}(Int(m), Int(n), colptr, rowval, label, nzval)
        check && check_invariants(d)
        return d
    end
end

LabelledSparseMatrix(m::Integer, n::Integer, colptr::Vector{Int}, rowval::Vector{Int32},
                     label::Vector{L}, nzval::Vector{Tv}; check::Bool = true) where {L, Tv} =
    LabelledSparseMatrix{L, Tv}(m, n, colptr, rowval, label, nzval; check)

function LabelledSparseMatrix{L, Tv}(m::Integer, n::Integer) where {L, Tv}
    m >= 0 && n >= 0 || throw(ArgumentError("negative size $m×$n"))
    return LabelledSparseMatrix{L, Tv}(m, n, ones(Int, Int(n) + 1), Int32[], L[], Tv[])
end

_label_valid(l::Integer) = l >= 0
_label_valid(l::NTuple{2, Integer}) = l[1] >= 0 && l[2] >= 0

"""
    check_invariants(d::LabelledSparseMatrix)

Check the invariants of `d` in O(nnz) (see [`LabelledSparseMatrix`](@ref)); throw an
`ArgumentError` naming the first violation, return `d` otherwise.
"""
function check_invariants(d::LabelledSparseMatrix)
    m, n = d.m, d.n
    m >= 0 && n >= 0 || throw(ArgumentError("negative size $m×$n"))
    m <= typemax(Int32) || throw(ArgumentError("$m rows: at most typemax(Int32) allowed"))
    length(d.colptr) == n + 1 ||
        throw(ArgumentError("colptr has length $(length(d.colptr)), expected $(n + 1)"))
    k = length(d.rowval)
    length(d.label) == k && length(d.nzval) == k ||
        throw(ArgumentError("rowval, label and nzval have lengths $k, $(length(d.label)), $(length(d.nzval))"))
    d.colptr[1] == 1 || throw(ArgumentError("colptr[1] = $(d.colptr[1]), expected 1"))
    d.colptr[n + 1] == k + 1 ||
        throw(ArgumentError("colptr[end] = $(d.colptr[n + 1]), expected $(k + 1)"))
    for j in 1:n
        a, b = d.colptr[j], d.colptr[j + 1] - 1
        a <= b + 1 && b <= k || throw(ArgumentError("colptr is not nondecreasing at column $j"))
        for p in a:b
            r = d.rowval[p]
            1 <= r <= m || throw(ArgumentError("row $r out of 1:$m (column $j, entry $p)"))
            _label_valid(d.label[p]) ||
                throw(ArgumentError("invalid label $(d.label[p]) (column $j, entry $p)"))
            iszero(d.nzval[p]) && throw(ArgumentError("stored zero (column $j, entry $p)"))
            if p > a
                r0 = d.rowval[p - 1]
                r0 < r || (r0 == r && isless(d.label[p - 1], d.label[p])) ||
                    throw(ArgumentError("(row, label) not strictly increasing (column $j, entry $p)"))
            end
        end
    end
    return d
end

"""
    labelled_sparse(I, J, L, V, m, n)

The `m × n` [`LabelledSparseMatrix`](@ref) whose monomials are `V[k] · L[k]` in row
`I[k]` and column `J[k]`, in any order (arguments in the order of `SparseArrays.sparse`,
labels before values). Monomials with the same row, column and label are summed
(checked arithmetic for machine integers), zero sums are dropped. Cost: a counting sort
on the columns, O(k + n), then a sort inside each column.
"""
function labelled_sparse(I::AbstractVector{<:Integer}, J::AbstractVector{<:Integer},
                         L::AbstractVector{LT}, V::AbstractVector{Tv},
                         m::Integer, n::Integer) where {LT, Tv}
    k = length(I)
    length(J) == k && length(L) == k && length(V) == k ||
        throw(ArgumentError("I, J, L and V have lengths $k, $(length(J)), $(length(L)), $(length(V))"))
    m >= 0 && n >= 0 || throw(ArgumentError("negative size $m×$n"))
    m <= typemax(Int32) || throw(ArgumentError("$m rows: at most typemax(Int32) allowed"))

    # Column pointers of the unmerged monomials (counting sort on the columns).
    ptr = zeros(Int, n + 1)
    for p in 1:k
        i, j = I[p], J[p]
        1 <= j <= n || throw(ArgumentError("column $j out of 1:$n (monomial $p)"))
        1 <= i <= m || throw(ArgumentError("row $i out of 1:$m (monomial $p)"))
        _label_valid(L[p]) || throw(ArgumentError("invalid label $(L[p]) (monomial $p)"))
        ptr[j + 1] += 1
    end
    ptr[1] = 1
    for j in 1:n
        ptr[j + 1] += ptr[j]
    end
    rows = Vector{Int32}(undef, k)
    labs = Vector{LT}(undef, k)
    vals = Vector{Tv}(undef, k)
    next = ptr[1:n]
    for p in 1:k
        j = J[p]
        q = next[j]
        next[j] = q + 1
        rows[q] = I[p]
        labs[q] = L[p]
        vals[q] = V[p]
    end

    # Sort each column by (row, label), merge equal keys, compact in place.
    colptr = Vector{Int}(undef, n + 1)
    colptr[1] = 1
    buf = Tuple{Int32, LT, Tv}[]
    w = 0
    for j in 1:n
        empty!(buf)
        for p in ptr[j]:(ptr[j + 1] - 1)
            push!(buf, (rows[p], labs[p], vals[p]))
        end
        sort!(buf; by = t -> (t[1], t[2]))
        q = 1
        while q <= length(buf)
            r, l, v = buf[q]
            q += 1
            while q <= length(buf) && buf[q][1] == r && buf[q][2] == l
                v = _checked_add(v, buf[q][3])
                q += 1
            end
            if !iszero(v)
                w += 1
                rows[w] = r
                labs[w] = l
                vals[w] = v
            end
        end
        colptr[j + 1] = w + 1
    end
    resize!(rows, w)
    resize!(labs, w)
    resize!(vals, w)
    return LabelledSparseMatrix{LT, Tv}(m, n, colptr, rows, labs, vals; check = false)
end

AbstractAlgebra.number_of_rows(d::LabelledSparseMatrix) = d.m
AbstractAlgebra.number_of_columns(d::LabelledSparseMatrix) = d.n
Base.size(d::LabelledSparseMatrix) = (d.m, d.n)
Base.size(d::LabelledSparseMatrix, i::Integer) = i == 1 ? d.m : i == 2 ? d.n : 1
SparseArrays.nnz(d::LabelledSparseMatrix) = length(d.nzval)
SparseArrays.nzrange(d::LabelledSparseMatrix, j::Integer) = d.colptr[j]:(d.colptr[j + 1] - 1)
SparseArrays.rowvals(d::LabelledSparseMatrix) = d.rowval
SparseArrays.nonzeros(d::LabelledSparseMatrix) = d.nzval

"""
    labels(d::LabelledSparseMatrix)

The labels of the stored entries of `d`, in storage order (as `rowvals` and `nonzeros`).
"""
labels(d::LabelledSparseMatrix) = d.label

Base.iszero(d::LabelledSparseMatrix) = isempty(d.nzval)

Base.:(==)(d::LabelledSparseMatrix, e::LabelledSparseMatrix) =
    d.m == e.m && d.n == e.n && d.colptr == e.colptr && d.rowval == e.rowval &&
    d.label == e.label && d.nzval == e.nzval
Base.hash(d::LabelledSparseMatrix, h::UInt) =
    hash(d.nzval, hash(d.label, hash(d.rowval, hash(d.colptr, hash((d.m, d.n), h)))))

function Base.show(io::IO, d::LabelledSparseMatrix{L, Tv}) where {L, Tv}
    k = nnz(d)
    print(io, d.m, "×", d.n, " LabelledSparseMatrix{", L, ", ", Tv, "} with ", k,
          k == 1 ? " stored entry" : " stored entries")
end

"""
    specialise(d::LabelledSparseMatrix)
    specialise(d::LabelledSparseMatrix, χ)

Integer matrix obtained by sending each label `γ` to the integer `χ(γ)`, or to `1` if
`χ` is omitted: the monomials of a given row and column are then summed, whatever
their labels (trivial module on both sides). One pass over the stored entries; the
result is a `SparseMatrixCSC{Tv, Int}` without stored zeros, computed with checked
arithmetic for machine integers. `χ` is any callable, for instance `(σ, τ) ↦ (−1)^ℓ(τ)`,
the sign module on the right (ℤε), with lengths given by the germ. The sign on both sides,
`(σ, τ) ↦ (−1)^(ℓ(σ) + ℓ(τ))`, only changes the signs of the cells of a homogeneous
resolution (`ℓ(σ) + ℓ(τ)` is the difference of the lengths of the weights of the source and
of the target), hence gives the trivial homology. For a resolution, the typed interface is
`homology(R, M)` and `Gauss.tensor_product`, with a module of coefficients `M`.
"""
specialise(d::LabelledSparseMatrix) = _specialise(d, nothing)
specialise(d::LabelledSparseMatrix, χ) = _specialise(d, χ)

_weight(::Nothing, l, v) = v
_weight(χ, l, v::Tv) where {Tv} = _checked_mul(convert(Tv, χ(l)), v)

function _specialise(d::LabelledSparseMatrix{L, Tv}, χ) where {L, Tv}
    colptr = Vector{Int}(undef, d.n + 1)
    colptr[1] = 1
    rowval = Int[]
    nzval = Tv[]
    sizehint!(rowval, nnz(d))
    sizehint!(nzval, nnz(d))
    for j in 1:d.n
        p, b = d.colptr[j], d.colptr[j + 1] - 1
        while p <= b
            r = d.rowval[p]
            acc = _weight(χ, d.label[p], d.nzval[p])
            p += 1
            while p <= b && d.rowval[p] == r
                acc = _checked_add(acc, _weight(χ, d.label[p], d.nzval[p]))
                p += 1
            end
            if !iszero(acc)
                push!(rowval, r)
                push!(nzval, acc)
            end
        end
        colptr[j + 1] = length(rowval) + 1
    end
    return SparseMatrixCSC{Tv, Int}(d.m, d.n, colptr, rowval, nzval)
end
