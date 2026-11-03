# The detailed displays (J1.h, point 5b): a complex shows its bases, one line per degree, and
# the boundaries of its generators; its differential the boundaries alone; an element, in
# full. With :limit (the REPL), everything keeps within the screen: the lines of the bases
# are cut in the middle, the list of the boundaries keeps its beginning and its end, a long
# boundary is cut at the width of the screen; for bases of cells, the boundaries of degree 1,
# d(a) = a[] - []a, are omitted but for a reminder. Without :limit, everything is written.

# Hooks of the bases: the title of the section of the bases, whether the boundaries of the
# lowest degree but one are omitted (cells), and the writing of a generator without brackets.
_basis_title(::AbstractVector) = "Basis"  # COV_EXCL_LINE (a constant, folded by the compiler)
_omits_first_boundaries(::AbstractVector) = false  # COV_EXCL_LINE (a constant, folded by the compiler)
_write_plain_generator(io::IO, ::Notation, b::AbstractVector, i::Integer) = print(io, b[i])

# Writes v, cut after the last term that keeps it within `width` characters, followed by
# " … (k terms)" (the first term is always written).
function _write_element_limited(io::IO, v::BasedFreeModuleElement, w::_ElementWriter, width::Integer)
    iszero(v) && return print(io, "0")
    n = length(v.indices)
    tail = " … ($n terms)"
    used = 0
    for k in 1:n
        piece = sprint(context = io) do s
            c = v.coefficients[k]
            print(s, k == 1 ? (c < 0 ? "-" : "") : (c < 0 ? " - " : " + "))
            _write_term(s, w, v.parent, c, v.indices[k], v.labels[k])
        end
        reserve = k < n ? textwidth(tail) : 0
        k > 1 && used + textwidth(piece) + reserve > width && return print(io, tail)
        print(io, piece)
        used += textwidth(piece)
    end
    return nothing
end

# With :limit, on one line (inside a vector, for instance); at the REPL, within the screen.
function Base.show(io::IO, v::BasedFreeModuleElement)
    w = _ElementWriter(v.parent, notation(io))
    get(io, :limit, false)::Bool ? _write_element_limited(io, v, w, displaysize(io)[2]) : _write_element(io, v, w)
    return nothing
end

# Inside a container, which sets :typeinfo (the arrays print their elements with this method,
# d.(basis(P[n]))), on one line; alone, within the screen.
function Base.show(io::IO, ::MIME"text/plain", v::BasedFreeModuleElement)
    w = _ElementWriter(v.parent, notation(io))
    get(io, :limit, false)::Bool || return _write_element(io, v, w)
    rows, columns = displaysize(io)
    width = get(io, :typeinfo, nothing) === nothing ? max(rows - 4, 1) * columns : columns - 2
    _write_element_limited(io, v, w, width)
    return nothing
end

# What the lines of the boundaries of the degree n need, computed once.
struct _BoundaryWriter{G, W <: _ElementWriter}
    generators::G             # for the generators of C_n, without brackets
    elements::W               # for the elements of C_{n-1}
end

_BoundaryWriter(C::AbstractChainComplex, N::Notation, n::Integer) =
    _BoundaryWriter(_generator_writer(C[n].basis, N), _ElementWriter(C[n - 1], N))

# The line "d(x) = …" of the i-th generator of degree n, its boundary within `width`
# characters (typemax(Int): in full).
function _write_boundary(io::IO, C::AbstractChainComplex, N::Notation, w::_BoundaryWriter, n::Integer,
                         i::Integer, width::Integer)
    head = sprint(context = io) do s
        print(s, N.differential_symbol, "(")
        _write_plain_generator(s, w.generators, C[n].basis, i)
        print(s, ") = ")
    end
    print(io, head)
    v = _apply(C, n, _element_of_generator(C[n], i))
    width == typemax(Int) ? _write_element(io, v, w.elements) :
                            _write_element_limited(io, v, w.elements, width - textwidth(head))
    return nothing
end

# The lines of the boundaries of C, `used` lines being already written.
function _show_boundaries(io::IO, C::AbstractChainComplex, used::Integer)
    N = notation(io)
    limited = get(io, :limit, false)::Bool
    rows, columns = displaysize(io)
    width = limited ? columns - 2 : typemax(Int)
    lo, hi = first(range(C)), last(range(C))
    writers = Dict(n => _BoundaryWriter(C, N, n) for n in (lo + 1):hi)
    start = limited && _omits_first_boundaries(C[lo].basis) ? lo + 2 : lo + 1
    lines = used
    if start > lo + 1 && rank(C[lo + 1]) > 0
        print(io, "\n  (boundaries of degree $(lo + 1) omitted, such as ")
        _write_boundary(io, C, N, writers[lo + 1], lo + 1, 1, typemax(Int))
        print(io, ")")
        lines += 1
    end
    degrees = start:hi
    counts = [rank(C[n]) for n in degrees]
    total = sum(counts; init = 0)
    budget = limited ? max(rows - lines - 2, 3) : total
    if total <= budget
        for k in 1:total
            _boundary_line(io, C, N, writers, degrees, counts, width, k)
        end
    else
        head = cld(budget - 1, 2)
        tail = budget - 1 - head
        for k in 1:head
            _boundary_line(io, C, N, writers, degrees, counts, width, k)
        end
        print(io, "\n  ⋮ ($(total - head - tail) more boundaries)")
        for k in (total - tail + 1):total
            _boundary_line(io, C, N, writers, degrees, counts, width, k)
        end
    end
    return nothing
end

# The line of the k-th boundary of the degrees `degrees`, of ranks `counts`.
function _boundary_line(io::IO, C::AbstractChainComplex, N::Notation, writers, degrees, counts::Vector{Int},
                        width::Integer, k::Integer)
    n, i = _locate(degrees, counts, k)
    print(io, "\n  ")
    _write_boundary(io, C, N, writers[n], n, i, width)
    return nothing
end

# The degree and the index of the k-th generator of the degrees `degrees`, of ranks `counts`.
function _locate(degrees, counts::Vector{Int}, k::Integer)
    r = Int(k)
    for (n, c) in zip(degrees, counts)
        r <= c && return (n, r)
        r -= c
    end
    throw(BoundsError(counts, k))    # COV_EXCL_LINE (k is at most the sum of the counts)
end

function Base.show(io::IO, ::MIME"text/plain", C::AbstractChainComplex)
    show(io, C)
    print(io, "\n", _basis_title(C[first(range(C))].basis), ":")
    width = displaysize(io)[2]
    for n in range(C)
        prefix = "  $n: "
        print(io, "\n", prefix)
        _print_limited(io, basis(C[n]), width - textwidth(prefix))
    end
    if length(range(C)) > 1
        print(io, "\nDifferential:")
        _show_boundaries(io, C, length(range(C)) + 3)
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", D::Differential)
    show(io, D)
    length(range(D.complex)) > 1 && _show_boundaries(io, D.complex, 1)
    return nothing
end
