# The elements of a germ for the user (J1.h): a simple with its germ, written as a word in
# the atoms, built from a word or read from its written form. The construction keeps its
# numbers (Int32); only what is handed to the user is wrapped.

"""
    GermElement

An element of a germ `G` (a simple), for the user: the germ and the number of the simple.
It is written as a word in the atoms, following the settings of the display
([`Notation`](@ref)): `121` for B⁺(A₂), `e` for the identity, `1·12·3` in a germ with at
least 10 atoms.

Built by calling the germ:
- `G()`, the identity; `G(a)`, the atom `a`; `G(1, 2, 1)`, the product of the atoms of a
  word, an `ArgumentError` if it is not a simple;
- `G("121")`, read from its written form, with the settings of the session; spaces, `·`,
  `⋅` and `.` between the atoms are accepted, and parentheses around the word.

Also given by [`simples`](@ref)`(G)` and [`elements`](@ref)`(G)`, and by the accessors of
the cells and of the elements of free modules. Operations: [`word`](@ref)`(s)`, the
smallest word in the atoms; `length(s)`; `isone(s)`, `one(s)`; `s * t`, the product if it
is a simple (an `ArgumentError` otherwise: the products of the monoid, by normal forms, are
not supported yet); `s \\ t`, the left quotient if `s ≼ t`; equality, hashing and the
degree-lexicographic order (`isless`). A `GermElement` is a scalar for broadcasting.

`GermElement(G, s)` wraps the number `s` of a simple of `G` (checked against the simples
numbered by `G`): a tool of the developer.
"""
struct GermElement{G <: AbstractGerm}
    germ::G
    index::Int32
    function GermElement(G::AbstractGerm, s::Integer; check::Bool = true)
        check && _check_simple(G, s)
        return new{typeof(G)}(G, Int32(s))
    end
end

# The number of simples numbered by G: all of them for a tabulated germ, those met so far
# for a lazy germ that counts them, nothing if unknown.
_numbered_count(G::AbstractGerm) =
    applicable(_simple_numbers, G) ? length(_simple_numbers(G)) :
    applicable(nnumbered, G) ? Int(nnumbered(G)) : nothing

function _check_simple(G::AbstractGerm, s::Integer)
    n = _numbered_count(G)
    0 <= s && (n === nothing || s < n) ||
        throw(ArgumentError("GermElement: $s is not the number of a simple of the germ" *
                            (n === nothing ? "" : " (0:$(n - 1))")))
    return nothing
end

# ------------------------------------------------------------------------------------
# Construction by calling the germ

(G::AbstractGerm)() = GermElement(G, 0; check = false)

function (G::AbstractGerm)(a::Integer)
    N = natoms(G)
    1 <= a <= N && return GermElement(G, a; check = false)
    throw(ArgumentError(_atom_message(N, a)))
end

function (G::AbstractGerm)(a::Integer, b::Integer, rest::Integer...)
    w = (a, b, rest...)
    return GermElement(G, _word_product(G, w, "G(" * join(w, ", ") * ")"); check = false)
end

(G::AbstractGerm)(text::AbstractString) = GermElement(G, _parse_simple(G, notation(), text); check = false)

# The message for G(a) with a not an atom: G(121) is the atom 121, not the word 1·2·1.
function _atom_message(N::Int, a::Integer)
    d = a > 0 ? reverse(digits(a)) : Int[]
    msg = "G($a): no atom $a, the atoms of the germ are 1:$N"
    if length(d) >= 2 && all(x -> 1 <= x <= N, d)
        written = N <= 9 ? join(d) : join(d, "·")
        msg *= "; for the word $(join(d, '·')), write G($(join(d, ", "))) or G(\"$written\")"
    end
    return msg
end

# The product of the atoms of the word w, which must be a simple; `context` names the call
# in the messages.
function _word_product(G::AbstractGerm, w, context::AbstractString)
    N = natoms(G)
    s = Int32(0)
    for a in w
        1 <= a <= N || throw(ArgumentError("$context: no atom $a, the atoms of the germ are 1:$N"))
        p = iszero(s) ? Int32(a) : product(G, s, a)
        p === nothing && throw(ArgumentError("$context: the word $(join(w, '·')) is not a simple of the germ"))
        s = Int32(p)
    end
    return s
end

# ------------------------------------------------------------------------------------
# Operations

function _same_germ(s::GermElement, t::GermElement, op::AbstractString)
    s.germ === t.germ || throw(ArgumentError("$op: the two elements belong to different germs"))
    return nothing
end

"""
    word(s::GermElement) -> Vector{Int32}

The smallest reduced word of `s` in the atoms of its germ, for the lexicographic order:
`word(s.germ, s.index)`. Empty for the identity.
"""
word(s::GermElement) = word(s.germ, s.index)

Base.length(s::GermElement) = length(s.germ, s.index)
Base.isone(s::GermElement) = iszero(s.index)
Base.one(s::GermElement) = GermElement(s.germ, 0; check = false)
Base.broadcastable(s::GermElement) = Ref(s)

function Base.:*(s::GermElement, t::GermElement)
    _same_germ(s, t, "*")
    G = s.germ
    p = iszero(s.index) ? t.index : iszero(t.index) ? s.index : product(G, s.index, t.index)
    p === nothing && throw(ArgumentError("$s * $t is not a simple of the germ: the products of elements " *
                                         "of the monoid (normal forms) are not supported yet"))
    return GermElement(G, p; check = false)
end

function Base.:\(s::GermElement, t::GermElement)
    _same_germ(s, t, "\\")
    q = ldiv(s.germ, s.index, t.index)
    q === nothing && throw(ArgumentError("$s \\ $t: $s does not left-divide $t"))
    return GermElement(s.germ, q; check = false)
end

Base.:(==)(s::GermElement, t::GermElement) = s.germ === t.germ && s.index == t.index
Base.hash(s::GermElement, h::UInt) = hash(s.index, hash(objectid(s.germ), hash(:GermElement, h)))

function Base.isless(s::GermElement, t::GermElement)
    _same_germ(s, t, "isless")
    return isless_deglex(s.germ, s.index, t.index)
end

# ------------------------------------------------------------------------------------
# The simples and the elements of a germ, as GermElements

# The simples of a germ whose simples are all numbered, created on demand.
struct _GermElements{G <: AbstractGerm} <: AbstractVector{GermElement{G}}
    germ::G
    n::Int
end

Base.size(v::_GermElements) = (v.n,)
Base.IndexStyle(::Type{<:_GermElements}) = IndexLinear()
Base.@propagate_inbounds function Base.getindex(v::_GermElements, i::Int)
    @boundscheck checkbounds(v, i)
    return GermElement(v.germ, i - 1; check = false)
end
Base.showarg(io::IO, v::_GermElements, toplevel::Bool) = print(io, "simples of ", v.germ)

"""
    simples(G::AbstractGerm)

The simples of `G`, as [`GermElement`](@ref)s, in the numbering of `G` (the canonical
order for a tabulated germ: by length, then smallest word, so that the last one is Δ).
Defined for a germ whose simples are all numbered (a tabulated germ); an `ArgumentError`
for a lazy germ, which does not list its simples. The same as [`elements`](@ref)`(G)`.
"""
function simples(G::AbstractGerm)
    applicable(_simple_numbers, G) ||
        throw(ArgumentError("simples: a $(nameof(typeof(G))) does not list its simples (lazy germ)"))
    return _GermElements(G, length(_simple_numbers(G)))
end

"""
    elements(G::AbstractGerm)

The elements of the germ `G`, as [`GermElement`](@ref)s: the same as
[`simples`](@ref)`(G)`, under the name of the germs, which need not come from a Garside
structure. Public, not exported (PermGroups, hence Chevie, and FiniteFields, Hecke and Oscar
export an `elements` of their own): with Chevie loaded, the `elements` of Chevie has a
method for germs (package extension), and `elements(G)` works without qualification.
"""
elements(G::AbstractGerm) = simples(G)

# ------------------------------------------------------------------------------------
# The Garside element

# The Garside element of G, the maximal simple: the last simple of a tabulated germ; for
# another germ, reached from the identity by right multiplications by atoms (each simple
# other than Δ has an atom a with σa simple, since σ \ Δ has a left descent).
function _garside_element(G::AbstractGerm)
    N = natoms(G)
    s = Int32(0)
    while true
        a = 1
        while a <= N
            p = iszero(s) ? Int32(a) : product(G, s, a)
            if p !== nothing
                s = Int32(p)
                break
            end
            a += 1
        end
        a > N && return s
    end
end

_garside_element(G::TableGerm) = Int32(nsimples(G) - 1)

# ------------------------------------------------------------------------------------
# Writing

# What the writing of the simples of a germ needs, computed once for a display: the names
# of the atoms, their separator, the identity, and the Garside element with its symbol (-1
# when it is written as a word).
struct _SimpleWriter{G <: AbstractGerm}
    germ::G
    names::Vector{String}
    separator::String
    identity::String
    garside::Int32
    garside_symbol::String
end

function _SimpleWriter(G::AbstractGerm, N::Notation)
    names = Complexes._atom_names(N, natoms(G))
    Δ = N.garside_symbol === nothing ? Int32(-1) : _garside_element(G)
    return _SimpleWriter(G, names, Complexes._atom_separator(N, names), N.identity, Δ,
                         something(N.garside_symbol, ""))
end

function _write_simple(io::IO, w::_SimpleWriter, s::Integer)
    iszero(s) && return print(io, w.identity)
    s == w.garside && return print(io, w.garside_symbol)
    for (i, a) in enumerate(word(w.germ, s))
        i > 1 && print(io, w.separator)
        print(io, w.names[a])
    end
    return nothing
end

Base.show(io::IO, s::GermElement) = _write_simple(io, _SimpleWriter(s.germ, notation(io)), s.index)

# ------------------------------------------------------------------------------------
# Reading

# Accepted between the atoms of a word, whatever the separator of the display.
const _ATOM_SEPARATORS = ('·', '⋅', '.', '•', ' ', '\t')

# The text without the parentheses that enclose all of it.
function _strip_parentheses(text::AbstractString)
    t = strip(text)
    while startswith(t, '(') && endswith(t, ')') && _closes_at_end(t)
        t = strip(t[nextind(t, firstindex(t)):prevind(t, lastindex(t))])
    end
    return t
end

# Whether the parenthesis that opens t closes at its end, as in (12), not in (1)(2).
function _closes_at_end(t::AbstractString)
    depth = 0
    for k in eachindex(t)
        c = t[k]
        depth += c == '(' ? 1 : c == ')' ? -1 : 0
        depth == 0 && return k == lastindex(t)
    end
    return false
end

# Splits the text before each character that is not a subscript digit or a prime: the
# names of a germ whose atoms are written side by side.
function _split_names(t::AbstractString)
    tokens = String[]
    buffer = IOBuffer()
    for c in t
        if !Complexes._is_modifier(c) && position(buffer) > 0
            push!(tokens, String(take!(buffer)))
        end
        print(buffer, c)
    end
    position(buffer) > 0 && push!(tokens, String(take!(buffer)))
    return tokens
end

# The number of the simple written `text` in the germ G with the settings N; `context` is
# the expression read, for the messages.
function _parse_simple(G::AbstractGerm, N::Notation, text::AbstractString; context::AbstractString = text)
    t = _strip_parentheses(text)
    isempty(t) && throw(ArgumentError("cannot read \"$context\": an empty simple"))
    t == N.identity && return Int32(0)
    N.garside_symbol !== nothing && t == N.garside_symbol && return _garside_element(G)
    names = Complexes._atom_names(N, natoms(G))
    if N.atom_separator isa String && !isempty(N.atom_separator)
        t = replace(t, N.atom_separator => " ")
    end
    tokens = any(in(_ATOM_SEPARATORS), t) ? split(t, collect(_ATOM_SEPARATORS); keepempty = false) :
             all(Complexes._is_self_delimiting, names) ? _split_names(t) : [t]
    index = Dict{String, Int}(names[a] => a for a in eachindex(names))
    w = Int[_atom_of(index, names, tok, context) for tok in tokens]
    return _word_product(G, w, "cannot read \"$context\"")
end

# The atom named `tok`, or numbered `tok`.
function _atom_of(index::Dict{String, Int}, names::Vector{String}, tok::AbstractString, context::AbstractString)
    a = get(index, tok, 0)
    a > 0 && return a
    n = tryparse(Int, tok)
    n !== nothing && 1 <= n <= length(names) && return n
    throw(ArgumentError("cannot read \"$context\": no atom named \"$tok\" (the germ has $(length(names)) atoms)"))
end

# ------------------------------------------------------------------------------------
# The germs as monoids of the modules of the complexes layer (J1.h): their elements are the
# labels of the terms, written as words, read from their written form, multiplied in the
# germ (a product that is not a simple is not a label), and handed to the user as
# GermElements, which act on the sides of the elements of the modules.

Complexes._label_writer(G::AbstractGerm, N::Notation) = _SimpleWriter(G, N)
Complexes._write_label(io::IO, w::_SimpleWriter, l::Integer) = _write_simple(io, w, l)
Complexes._digit_labels(w::_SimpleWriter) = any(name -> isdigit(first(name)), w.names)
Complexes._read_label(G::AbstractGerm, N::Notation, text::AbstractString, context::AbstractString) =
    _parse_simple(G, N, text; context = context)
Complexes._label_product(G::AbstractGerm, a::Integer, b::Integer) =
    iszero(a) ? Int32(b) : iszero(b) ? Int32(a) : product(G, a, b)
Complexes._label_value(G::AbstractGerm, l::Int32) = GermElement(G, l; check = false)

function Complexes._label_code(G::AbstractGerm, s::GermElement)
    s.germ === G || throw(ArgumentError("the element $s is not one of the germ of the module"))
    return s.index
end

Base.:*(s::GermElement, v::Complexes.BasedFreeModuleElement) = Complexes._act(v, :left, s.germ, s.index)
Base.:*(v::Complexes.BasedFreeModuleElement, s::GermElement) = Complexes._act(v, :right, s.germ, s.index)
