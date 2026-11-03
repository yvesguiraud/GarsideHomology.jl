# Elements of the free modules with a basis (J1.h): sparse vectors over the basis, with the
# labels of the monoid algebra, written as sums of terms c×σ[x]τ and read from that form;
# the differential of a complex as an object applied to them.

# ------------------------------------------------------------------------------------
# Hooks of the monoids: writing, reading, multiplying and handing over the labels. The
# defaults number the elements (0 is the identity) and multiply nothing but the identity;
# the germs of GarsideHomology write them as words in the atoms and multiply them in the
# germ.

# The writer of the labels of a monoid with the settings N, computed once for a display.
struct _NumberedLabels
    identity::String
end
_label_writer(::AbstractMonoid, N::Notation) = _NumberedLabels(N.identity)
_write_label(io::IO, w::_NumberedLabels, l::Integer) = iszero(l) ? print(io, w.identity) : print(io, l)

# Whether the written labels can begin with a digit, so that a coefficient needs a marker.
_digit_labels(::_NumberedLabels) = true  # COV_EXCL_LINE (a constant, folded by the compiler)

# The label written `text` (`context`, the expression read, for the messages).
function _read_label(::AbstractMonoid, N::Notation, text::AbstractString, context::AbstractString)
    t = strip(text)
    t == N.identity && return Int32(0)
    l = tryparse(Int32, t)
    l !== nothing && l >= 0 && return l
    throw(ArgumentError("cannot read \"$context\": \"$t\" is not an element of the monoid"))
end

# The product of two labels, nothing when it is not a label (not a simple, for a germ).
_label_product(::AbstractMonoid, a::Integer, b::Integer) =
    iszero(a) ? Int32(b) : iszero(b) ? Int32(a) : nothing

# The label of an element of the monoid given by the user, and the element of a label.
_label_code(::AbstractMonoid, s::Integer) = Int32(s)
_label_value(::AbstractMonoid, l::Int32) = l

# ------------------------------------------------------------------------------------
# Hooks of the bases: writing a generator in brackets, finding it from itself or from its
# written form, and its degree. The defaults use the printed form of the descriptors.

_generator_writer(::AbstractVector, N::Notation) = N
_write_generator(io::IO, ::Notation, b::AbstractVector, i::Integer) = print(io, "[", b[i], "]")

# The index in the basis b of its element x.
function _basis_index(b::AbstractVector, x)
    i = findfirst(isequal(x), b)
    i === nothing && throw(ArgumentError("$(repr(x)) is not an element of the basis"))
    return Int(i)
end

# The index in the basis b of the element written `text`.
function _find_in_basis(b::AbstractVector, ::Notation, text::AbstractString, context::AbstractString)
    t = strip(text)
    i = findfirst(x -> string(x) == t, b)
    i === nothing && throw(ArgumentError("cannot read \"$context\": \"$t\" is not an element of the basis"))
    return Int(i)
end

# The degree that the basis b knows, nothing if it does not know one.
_basis_degree(::AbstractVector) = nothing  # COV_EXCL_LINE (a constant, folded by the compiler)

# The degree and the index of the generator of C written `text`, b0 being the basis of the
# first degree, by its type: by default, the element whose printed form is the text, in
# any degree.
function _find_generator(::AbstractVector, C::AbstractChainComplex, ::Notation, text::AbstractString)
    t = strip(text)
    found = [(n, i) for n in range(C) for (i, x) in enumerate(basis(C[n])) if string(x) == t]
    isempty(found) && throw(ArgumentError("\"$text\" is not a generator of the complex"))
    length(found) == 1 ||
        throw(ArgumentError("\"$text\" names generators in several places (degrees $(unique(first.(found))))"))
    return only(found)
end

# The degree of a generator x of C given by the user, b0 being the basis of the first
# degree: by default, a descriptor does not know its degree.
_generator_degree(::AbstractVector, x) =
    throw(ArgumentError("$(repr(x)) does not know its degree: apply the differential to C[n](x)"))

# ------------------------------------------------------------------------------------
# The elements

"""
    BasedFreeModuleElement

An element of a free module with a basis `M` ([`BasedFreeModule`](@ref)): a finite sum of
terms `c×σ[x]τ`, with `x` an element of the basis, `c` an integer and, over a monoid algebra,
`σ` and `τ` elements of the monoid on the sides of the module (both for bimodules, one for
left or right modules, none over a ring). Kept as sorted vectors of indices, labels and
coefficients (`Int`, with checked arithmetic), without zeros, so that equality is that of
the vectors. The terms are in the order of the basis, then, for bimodules, of `τ`, then of
`σ`: the terms `σ[x]` come before the terms `[x]τ`.

Built by calling the module: `M(x)`, the generator `[x]` for an element `x` of the basis (a
[`Cell`](@ref) for a resolution); `M("2×12[1|21]3 - [2|12]")`, read from its written form
(a text without brackets is a generator, `M("2|12")`); over a ring, `M(v)` from the vector of
its coefficients on the basis; `zero(M)`. For a complex, `C(x)` for a cell and
`C("…")` find the module themselves.

Operations: `+`, `-`, multiples by integers, the actions of the elements of the monoid on
the sides of the module (`σ * v`, `v * τ`, with [`GermElement`](@ref)s; an `ArgumentError`
if a product is not a simple), equality and hashing, `iszero`, `parent(v)` (the module);
`v[x]` (over a ring), `v[σ, x]`, `v[x, τ]` or `v[σ, x, τ]`, the coefficient of a term; the
iteration over the terms, as tuples `(c, σ, x, τ)` (`nothing` for an absent side);
`length(v)`, the number of terms. Over a ring, `Vector(v)` and `SparseVector(v)` give the
coefficients on the basis, of length `rank(M)`.

Written as `2×12[1|21]3 - [2|12]`, following the settings of the display
([`Notation`](@ref)): the identity is omitted on both sides, the coefficients ±1 too, the
cell of degree 0 is `[]`, and a marker `×` separates a coefficient from what follows when
the elements of the monoid are written with digits.
"""
struct BasedFreeModuleElement{M <: BasedFreeModule, L, Tv}
    parent::M
    indices::Vector{Int32}
    labels::Vector{L}
    coefficients::Vector{Tv}
end

# The type of the labels of the elements of M: nothing over a ring.
_element_label_type(::BasedFreeModule{<:AbstractAlgebra.Ring}) = Nothing  # COV_EXCL_LINE (a constant, folded by the compiler)
_element_label_type(M::BasedFreeModule{<:MonoidAlgebra}) = label_type(M.ring)

_identity_label(::Type{Nothing}) = nothing  # COV_EXCL_LINE (a constant, folded by the compiler)
_identity_label(::Type{L}) where {L} = identity_label(L)  # COV_EXCL_LINE (a constant, folded by the compiler)

# The order of the terms: by generator, then for bimodules by τ, then by σ, so that the
# terms σ[x] come before the terms [x]τ (d[a] = a[] - []a).
_order_key(i::Int32, ::Nothing) = (i, Int32(0), Int32(0))
_order_key(i::Int32, l::Int32) = (i, l, Int32(0))
_order_key(i::Int32, l::NTuple{2, Int32}) = (i, l[2], l[1])

# The element of M of the terms (I[k], L[k], V[k]), in canonical form: sorted by generator
# and label (see _order_key), equal terms summed (checked), zeros dropped.
function _element(M::BasedFreeModule, I::Vector{Int32}, L::Vector, V::Vector{Int})
    p = sortperm([_order_key(I[k], L[k]) for k in eachindex(I)])
    indices, labels, coefficients = Int32[], similar(L, 0), Int[]
    for k in p
        if !isempty(indices) && indices[end] == I[k] && labels[end] == L[k]
            coefficients[end] = Base.checked_add(coefficients[end], V[k])
        else
            push!(indices, I[k])
            push!(labels, L[k])
            push!(coefficients, V[k])
        end
    end
    keep = findall(!iszero, coefficients)
    return BasedFreeModuleElement(M, indices[keep], labels[keep], coefficients[keep])
end

_element_of_generator(M::BasedFreeModule, i::Integer) =
    BasedFreeModuleElement(M, Int32[i], [_identity_label(_element_label_type(M))], [1])

(M::BasedFreeModule)(x) = _element_of_generator(M, _basis_index(M.basis, x))
(M::BasedFreeModule)(text::AbstractString) = _read_element(M, text)

function (M::BasedFreeModule{<:AbstractAlgebra.Ring})(v::AbstractVector{<:Integer})
    length(v) == rank(M) ||
        throw(ArgumentError("a vector of $(length(v)) coefficients for a module of rank $(rank(M))"))
    keep = findall(!iszero, v)
    return BasedFreeModuleElement(M, Int32.(keep), fill(nothing, length(keep)), Int[v[k] for k in keep])
end

Base.zero(M::BasedFreeModule) = _element(M, Int32[], _element_label_type(M)[], Int[])
Base.zero(v::BasedFreeModuleElement) = zero(v.parent)
Base.iszero(v::BasedFreeModuleElement) = isempty(v.indices)
Base.parent(v::BasedFreeModuleElement) = v.parent
Base.length(v::BasedFreeModuleElement) = length(v.indices)

function _same_module(v::BasedFreeModuleElement, w::BasedFreeModuleElement, op::AbstractString)
    v.parent === w.parent || v.parent == w.parent ||
        throw(ArgumentError("$op: the two elements belong to different modules"))
    return nothing
end

function Base.:+(v::BasedFreeModuleElement, w::BasedFreeModuleElement)
    _same_module(v, w, "+")
    return _element(v.parent, vcat(v.indices, w.indices), vcat(v.labels, w.labels), vcat(v.coefficients, w.coefficients))
end

Base.:-(v::BasedFreeModuleElement) =
    BasedFreeModuleElement(v.parent, copy(v.indices), copy(v.labels), Base.checked_neg.(v.coefficients))
Base.:-(v::BasedFreeModuleElement, w::BasedFreeModuleElement) = (_same_module(v, w, "-"); v + (-w))

function Base.:*(c::Integer, v::BasedFreeModuleElement)
    iszero(c) && return zero(v)
    return BasedFreeModuleElement(v.parent, copy(v.indices), copy(v.labels), [Base.checked_mul(Int(c), a) for a in v.coefficients])
end
Base.:*(v::BasedFreeModuleElement, c::Integer) = c * v

Base.:(==)(v::BasedFreeModuleElement, w::BasedFreeModuleElement) =
    (v.parent === w.parent || v.parent == w.parent) && v.indices == w.indices && v.labels == w.labels &&
    v.coefficients == w.coefficients
Base.hash(v::BasedFreeModuleElement, h::UInt) =
    hash(v.coefficients, hash(v.labels, hash(v.indices, hash(rank(v.parent), hash(:BasedFreeModuleElement, h)))))

# The element σ·v (side :left) or v·τ (side :right) for the element s of the monoid.
function _act(v::BasedFreeModuleElement, side::Symbol, monoid, s::Int32)
    A = v.parent.ring
    A isa MonoidAlgebra || throw(ArgumentError("an element of a module over a ring has no action of the monoid"))
    A.monoid === monoid || throw(ArgumentError("the element of the monoid is not one of the monoid of the module"))
    A.side === :both || A.side === side ||
        throw(ArgumentError("a module of $(A.side) modules has no action on the $side"))
    labels = similar(v.labels)
    for k in eachindex(v.labels)
        labels[k] = _acted_label(A, side, v.labels[k], s)
    end
    return _element(v.parent, copy(v.indices), labels, copy(v.coefficients))
end

function _product_label(A::MonoidAlgebra, a::Integer, b::Integer)
    p = _label_product(A.monoid, a, b)
    p === nothing && throw(ArgumentError("a product of elements of the monoid is not a simple: the products of " *
                                         "the monoid (normal forms) are not supported yet"))
    return Int32(p)
end

_acted_label(A::MonoidAlgebra, side::Symbol, l::Int32, s::Int32) =
    side === :left ? _product_label(A, s, l) : _product_label(A, l, s)
_acted_label(A::MonoidAlgebra, side::Symbol, l::NTuple{2, Int32}, s::Int32) =
    side === :left ? (_product_label(A, s, l[1]), l[2]) : (l[1], _product_label(A, l[2], s))

# The coefficient of the term of index i and label l.
function _coefficient(v::BasedFreeModuleElement, i::Integer, l)
    for k in eachindex(v.indices)
        v.indices[k] == i && v.labels[k] == l && return v.coefficients[k]
    end
    return zero(eltype(v.coefficients))
end

Base.getindex(v::BasedFreeModuleElement{<:Any, Nothing}, x) = _coefficient(v, _basis_index(v.parent.basis, x), nothing)

function Base.getindex(v::BasedFreeModuleElement{<:Any, Int32}, a, b)
    A = v.parent.ring
    A.side === :left && return _coefficient(v, _basis_index(v.parent.basis, b), _label_code(A.monoid, a))
    return _coefficient(v, _basis_index(v.parent.basis, a), _label_code(A.monoid, b))
end

function Base.getindex(v::BasedFreeModuleElement{<:Any, NTuple{2, Int32}}, σ, x, τ)
    A = v.parent.ring
    return _coefficient(v, _basis_index(v.parent.basis, x), (_label_code(A.monoid, σ), _label_code(A.monoid, τ)))
end

# The sides of a label, for the user (nothing for an absent side).
_sides(::BasedFreeModule, ::Nothing) = (nothing, nothing)  # COV_EXCL_LINE (a constant, folded by the compiler)
function _sides(M::BasedFreeModule, l::Int32)
    s = _label_value(M.ring.monoid, l)
    return M.ring.side === :left ? (s, nothing) : (nothing, s)
end
_sides(M::BasedFreeModule, l::NTuple{2, Int32}) = (_label_value(M.ring.monoid, l[1]), _label_value(M.ring.monoid, l[2]))

function Base.iterate(v::BasedFreeModuleElement, k::Int = 1)
    k > length(v.indices) && return nothing
    σ, τ = _sides(v.parent, v.labels[k])
    return (v.coefficients[k], σ, v.parent.basis[v.indices[k]], τ), k + 1
end

function Base.Vector(v::BasedFreeModuleElement{<:Any, Nothing})
    out = zeros(eltype(v.coefficients), rank(v.parent))
    for k in eachindex(v.indices)
        out[v.indices[k]] = v.coefficients[k]
    end
    return out
end

SparseArrays.SparseVector(v::BasedFreeModuleElement{<:Any, Nothing}) =
    SparseArrays.sparsevec(Int.(v.indices), copy(v.coefficients), rank(v.parent))

# ------------------------------------------------------------------------------------
# Writing

# What the writing of the elements of a module needs, computed once for a display.
struct _ElementWriter{W, G}
    labels::W                # nothing over a ring
    generators::G
    marker::String           # after a coefficient other than ±1 ("" for none)
    parentheses::Bool        # the monomial in parentheses after such a coefficient
end

function _ElementWriter(M::BasedFreeModule, N::Notation)
    labels = M.ring isa MonoidAlgebra ? _label_writer(M.ring.monoid, N) : nothing
    left = M.ring isa MonoidAlgebra && M.ring.side !== :right
    style = N.coefficient_style
    marker = style === :times ? "×" : style === :star ? "*" :
             style === :auto && left && _digit_labels(labels) ? "×" : ""
    return _ElementWriter(labels, _generator_writer(M.basis, N), marker, style === :parentheses)
end

function _write_term(io::IO, w::_ElementWriter, M::BasedFreeModule, c::Integer, i::Integer, l)
    a = abs(c)
    if a != 1
        print(io, a)
        print(io, w.parentheses ? "(" : w.marker)
    end
    σ, τ = _label_parts(l, M)
    σ === nothing || iszero(σ) || _write_label(io, w.labels, σ)
    _write_generator(io, w.generators, M.basis, i)
    τ === nothing || iszero(τ) || _write_label(io, w.labels, τ)
    a != 1 && w.parentheses && print(io, ")")
    return nothing
end

_label_parts(::Nothing, ::BasedFreeModule) = (nothing, nothing)  # COV_EXCL_LINE (a constant, folded by the compiler)
_label_parts(l::Int32, M::BasedFreeModule) = M.ring.side === :left ? (l, nothing) : (nothing, l)
_label_parts(l::NTuple{2, Int32}, ::BasedFreeModule) = l

function _write_element(io::IO, v::BasedFreeModuleElement, w::_ElementWriter)
    iszero(v) && return print(io, "0")
    for k in eachindex(v.indices)
        c = v.coefficients[k]
        print(io, k == 1 ? (c < 0 ? "-" : "") : (c < 0 ? " - " : " + "))
        _write_term(io, w, v.parent, c, v.indices[k], v.labels[k])
    end
    return nothing
end

# The show methods are in display.jl, with the detailed displays.

# ------------------------------------------------------------------------------------
# Reading

# A term as written: its sign and coefficient, and the texts of σ, of the generator and of τ.
struct _TermText
    coefficient::Int
    left::String
    generator::String
    right::String
end

# The terms of `t` with their signs, cut at the signs + and - (or −) outside brackets and
# parentheses.
function _split_terms(t::AbstractString, context::AbstractString)
    terms = Tuple{Int, String}[]
    parentheses, brackets = 0, 0
    sign = 1
    buffer = IOBuffer()
    for c in t
        if parentheses == 0 && brackets == 0 && c in ('+', '-', '−')
            if position(buffer) > 0
                push!(terms, (sign, String(take!(buffer))))
                sign = 1
            end
            c == '+' || (sign = -sign)
            continue
        end
        parentheses += c == '(' ? 1 : c == ')' ? -1 : 0
        brackets += c == '[' ? 1 : c == ']' ? -1 : 0
        min(parentheses, brackets) < 0 &&
            throw(ArgumentError("cannot read \"$context\": a bracket or a parenthesis is not opened"))
        (isspace(c) && position(buffer) == 0) || print(buffer, c)
    end
    parentheses == 0 && brackets == 0 ||
        throw(ArgumentError("cannot read \"$context\": a bracket or a parenthesis is not closed"))
    position(buffer) > 0 || throw(ArgumentError("cannot read \"$context\": a term is missing"))
    push!(terms, (sign, String(take!(buffer))))
    return terms
end

# The index of the parenthesis that closes the one at the start of s.
function _closing_parenthesis(s::AbstractString)
    depth = 0
    for k in eachindex(s)
        depth += s[k] == '(' ? 1 : s[k] == ')' ? -1 : 0
        depth == 0 && return k
    end
    return nothing    # COV_EXCL_LINE (the terms have balanced parentheses, checked by _split_terms)
end

# The parts of a term: leading digits are a coefficient when a marker (×, *, ⋅) or a
# parenthesis follows them, or when the elements of the monoid on the left cannot begin with
# a digit (`digits` false); otherwise they begin σ, as in the writing.
function _term_text(body::AbstractString, digits::Bool, context::AbstractString)
    fail(msg) = throw(ArgumentError("cannot read \"$context\": $msg"))
    s = strip(body)
    coefficient = 1
    m = match(r"^[0-9]+", s)
    if m !== nothing
        rest = lstrip(s[(ncodeunits(m.match) + 1):end])
        if !isempty(rest) && first(rest) in ('×', '*', '⋅')
            coefficient = parse(Int, m.match)
            s = lstrip(rest[nextind(rest, firstindex(rest)):end])
        elseif startswith(rest, '(') || !digits
            coefficient = parse(Int, m.match)
            s = rest
        end
    end
    left = ""
    if startswith(s, '(')
        k = _closing_parenthesis(s)
        inner = strip(s[nextind(s, firstindex(s)):prevind(s, k)])
        after = strip(s[nextind(s, k):end])
        if occursin('[', inner)
            isempty(after) || fail("\"$after\" after a monomial in parentheses")
            s = inner
        else
            left = String(inner)
            startswith(after, '[') || fail("a generator in brackets must follow ($inner)")
            s = after
        end
    end
    k = findfirst('[', s)
    k === nothing && fail("\"$s\" has no generator in brackets")
    isempty(left) && (left = String(strip(s[firstindex(s):prevind(s, k)])))
    close = findnext(']', s, k)
    close === nothing && fail("a bracket is not closed")     # COV_EXCL_LINE (checked by _split_terms)
    generator = String(strip(s[nextind(s, k):prevind(s, close)]))
    right = String(strip(s[nextind(s, close):end]))
    return _TermText(coefficient, left, generator, right)
end

# Whether leading digits of a term over the ring R may begin σ rather than be a coefficient.
_digits_begin_labels(::AbstractAlgebra.Ring, ::Notation) = false  # COV_EXCL_LINE (a constant, folded by the compiler)
_digits_begin_labels(A::MonoidAlgebra, N::Notation) = A.side !== :right && _digit_labels(_label_writer(A.monoid, N))

# The label of a term over the ring R.
function _term_label(::AbstractAlgebra.Ring, ::Notation, t::_TermText, context::AbstractString)
    isempty(t.left) && isempty(t.right) ||
        throw(ArgumentError("cannot read \"$context\": a module over a ring has no element of a monoid on its sides"))
    return nothing
end

function _term_label(A::MonoidAlgebra, N::Notation, t::_TermText, context::AbstractString)
    σ = isempty(t.left) ? Int32(0) : _read_label(A.monoid, N, t.left, context)
    τ = isempty(t.right) ? Int32(0) : _read_label(A.monoid, N, t.right, context)
    A.side === :both && return (σ, τ)
    absent = A.side === :left ? t.right : t.left
    isempty(absent) ||
        throw(ArgumentError("cannot read \"$context\": a module of $(A.side) modules has no element of the monoid on the other side"))
    return A.side === :left ? σ : τ
end

# The terms of an expression with brackets, read with the settings N for the ring R.
function _read_terms(R, N::Notation, text::AbstractString)
    digits = _digits_begin_labels(R, N)
    return [(sign, _term_text(body, digits, text)) for (sign, body) in _split_terms(strip(text), text)]
end

function _read_element(M::BasedFreeModule, text::AbstractString)
    N = notation()
    t = strip(text)
    isempty(t) && throw(ArgumentError("cannot read an empty text"))
    occursin('[', t) || return t == "0" ? zero(M) : _element_of_generator(M, _find_in_basis(M.basis, N, t, text))
    L = _element_label_type(M)
    I, labels, V = Int32[], L[], Int[]
    for (sign, term) in _read_terms(M.ring, N, text)
        push!(I, _find_in_basis(M.basis, N, term.generator, text))
        push!(labels, _term_label(M.ring, N, term, text))
        push!(V, sign * term.coefficient)
    end
    return _element(M, I, labels, V)
end

"""
    (C::AbstractChainComplex)(text::AbstractString)

What `text` writes in the complex `C`: with brackets, an element of the module of the degree
of its generators (`P("2×12[1|21]3 - [2|12]")`, an `ArgumentError` if they are not all of the
same degree); without, the generator written `text`, in any degree: for a resolution of
GarsideHomology, the cell of bar form `text` (`P("2|12")`); for other bases, the element whose
printed form is `text` (an `ArgumentError` if none or several are, in several degrees).
"""
function (C::AbstractChainComplex)(text::AbstractString)
    N = notation()
    b0 = basis(C[first(range(C))])
    isempty(strip(text)) && throw(ArgumentError("cannot read an empty text"))
    if !occursin('[', text)
        n, i = _find_generator(b0, C, N, text)
        return basis(C[n])[i]
    end
    R = base_ring(C)
    terms = _read_terms(R, N, text)
    found = [_find_generator(b0, C, N, term.generator) for (_, term) in terms]
    degrees = unique(first.(found))
    length(degrees) == 1 ||
        throw(ArgumentError("cannot read \"$text\": its generators have the degrees $degrees"))
    M = C[only(degrees)]
    labels = [_term_label(R, N, term, text) for (_, term) in terms]
    return _element(M, Int32[i for (_, i) in found], labels, Int[sign * term.coefficient for (sign, term) in terms])
end

# An integer does not name a generator; C(2|12) is C(14) for Julia.
function (C::AbstractChainComplex)(n::Integer)
    throw(ArgumentError("C($n): an integer does not name a generator. If you wrote C(x|y), Julia " *
                        "computed the bitwise or x|y = $n before the call: write the bar form " *
                        "between quotes, C(\"x|y\")"))
end

# ------------------------------------------------------------------------------------
# The differential as an object

"""
    Differential

The differential `d` of a chain complex `C` as a whole, as returned by
[`differential`](@ref)`(C)`: `d(x)`, or `d[x]`, is the boundary of the generator `x` of
`C_{n+1}` (a [`Cell`](@ref) of a resolution, which knows its degree), and `d(v)` that of an
element `v` of a module of `C`, an element of `C_n`. For the terms `c×σ[x]τ` of a module
over a monoid algebra, `d(v)` is the sum of the `c × σ·d(x)·τ`, an `ArgumentError` if a
product of elements of the monoid is not a simple. The boundary of an element of the lowest
degree is 0, in a module of rank 0 below it (whose elements have the boundary 0 too).
`d(d(x))` checks `d² = 0` on `x`. At the REPL, `d` shows the boundaries of the generators,
within the screen; `d.(basis(C[n]))` gives those of one degree.
"""
struct Differential{C <: AbstractChainComplex}
    complex::C
end

"""
    differential(C::AbstractChainComplex) -> Differential

The differential of the complex `C` as a whole, a [`Differential`](@ref), to apply to its
generators and to the elements of its modules: `d = differential(P); d(P("2|12"))`.
"""
differential(C::AbstractChainComplex) = Differential(C)

Base.show(io::IO, D::Differential) = (print(io, "Differential of "); show(io, D.complex))

(D::Differential)(v::BasedFreeModuleElement) = _apply(D.complex, _degree_in(D.complex, v.parent), v)

function (D::Differential)(x)
    C = D.complex
    n = _generator_degree(basis(C[first(range(C))]), x)
    return D(C[n](x))
end

Base.getindex(D::Differential, x) = D(x)

# The degree of the module M in C.
function _degree_in(C::AbstractChainComplex, M::BasedFreeModule)
    d = _basis_degree(M.basis)
    candidates = d === nothing ? range(C) : (d in range(C) ? (d:d) : (1:0))
    found = [n for n in candidates if C[n] === M]
    isempty(found) && (found = [n for n in candidates if C[n] == M])
    # A module of rank 0 over the ring of C, not in C: that of the boundaries of the lowest
    # degree, below it, so that d(d(x)) is defined for every generator x.
    isempty(found) && rank(M) == 0 && M.ring == base_ring(C) && return first(range(C)) - 1
    isempty(found) && throw(ArgumentError("the element does not belong to a module of the complex"))
    length(found) == 1 || throw(ArgumentError("the module of the element is that of several degrees $found of the complex"))
    return only(found)
end

# The boundary of the element v of C_n, in C_{n-1}.
function _apply(C::AbstractChainComplex, n::Integer, v::BasedFreeModuleElement)
    target = (n - 1) in range(C) ? C[n - 1] : BasedFreeModule(base_ring(C), 0)
    I, L, V = Int32[], similar(v.labels, 0), Int[]
    _apply_terms!(I, L, V, differential(C, n - 1), v, base_ring(C))
    return _element(target, I, L, V)
end

function _apply_terms!(I, L, V, d::SparseMatrixCSC, v::BasedFreeModuleElement, ::AbstractAlgebra.Ring)
    rows, values = rowvals(d), nonzeros(d)
    for k in eachindex(v.indices)
        for p in nzrange(d, v.indices[k])
            push!(I, Int32(rows[p]))
            push!(L, nothing)
            push!(V, Base.checked_mul(v.coefficients[k], Int(values[p])))
        end
    end
    return nothing
end

function _apply_terms!(I, L, V, d::LabelledSparseMatrix, v::BasedFreeModuleElement, A::MonoidAlgebra)
    rows, labs, values = rowvals(d), labels(d), nonzeros(d)
    for k in eachindex(v.indices)
        l = v.labels[k]
        for p in nzrange(d, v.indices[k])
            push!(I, Int32(rows[p]))
            push!(L, _boundary_label(A, l, labs[p]))
            push!(V, Base.checked_mul(v.coefficients[k], Int(values[p])))
        end
    end
    return nothing
end

# The label of the term σ·(σ′[y]τ′)·τ: (σσ′, τ′τ) for bimodules, σσ′ on the left, τ′τ on
# the right.
_boundary_label(A::MonoidAlgebra, l::NTuple{2, Int32}, m::NTuple{2, Int32}) =
    (_product_label(A, l[1], m[1]), _product_label(A, m[2], l[2]))
_boundary_label(A::MonoidAlgebra, l::Int32, m::Int32) =
    A.side === :left ? _product_label(A, l, m) : _product_label(A, m, l)
