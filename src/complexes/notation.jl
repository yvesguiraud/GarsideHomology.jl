# Settings of the display (J1.h): how the elements of the monoids, the cells and the elements
# of the free modules are written, for the session or for one display.

"""
    Notation(; atom_names = nothing, atom_separator = :auto, identity = "e",
             garside_symbol = nothing, bar_separator = "|", empty_cell = "e",
             empty_cell_in_module = "[]", coefficient_style = :auto,
             differential_symbol = "d")
    Notation(N::Notation; settings...)

The settings of the display of the elements of a monoid (a simple is written as a word in
the atoms), of the cells (`x₁|⋯|xₙ`) and of the elements of the free modules (`c×σ[x]τ`).
The second form copies `N` with the settings given changed. The settings of the session
are changed by [`set_notation!`](@ref) and read by [`notation`](@ref); for one display,
give a `Notation` under the key `:notation` of an `IOContext`, which wins over the session:
`show(IOContext(stdout, :notation => Notation(atom_separator = ".")), x)`.

- `atom_names`: the names of the atoms `1:N`. `nothing`, their numbers; `:subscripts`,
  `a₁, a₂, …`; a vector of strings (the atoms beyond its length keep their numbers); or a
  function of the number of an atom.
- `atom_separator`: written between the atoms of a word. `:auto` writes them side by side
  when every name of an atom of the germ is a single character, possibly followed by
  subscript digits or primes (`1` to `9`, `a`, `s₁₂`, `b′`), and separates them by `·`
  otherwise, in the whole germ (germs with at least 10 atoms numbered); a string imposes
  the separator.
- `identity`: the identity of the monoid.
- `garside_symbol`: `nothing` writes the Garside element as a word; a string (`"Δ"`, `"δ"`)
  writes it by that symbol, for a germ that has one.
- `bar_separator`, `empty_cell`, `empty_cell_in_module`: the separator of the letters of a
  cell, the cell of degree 0, and its image `[e]` in a free module.
- `coefficient_style`: how an integer coefficient other than ±1 is separated from what
  follows it. `:auto` writes `×` only when the names of the atoms can begin with a digit,
  hence for the numbered atoms (`2×12[1|21]3`, `2×[x]`), and nothing otherwise (`2ab[x]c`);
  `:times` (`×`) and `:star` (`*`) always; `:parentheses` puts the monomial in parentheses
  (`2(12[1|21]3)`).
- `differential_symbol`: the name of the differential in the displays of complexes.

The reading of expressions follows the same settings.
"""
struct Notation
    atom_names::Union{Nothing, Symbol, Vector{String}, Function}
    atom_separator::Union{Symbol, String}
    identity::String
    garside_symbol::Union{Nothing, String}
    bar_separator::String
    empty_cell::String
    empty_cell_in_module::String
    coefficient_style::Symbol
    differential_symbol::String
    function Notation(atom_names, atom_separator, identity, garside_symbol, bar_separator,
                      empty_cell, empty_cell_in_module, coefficient_style, differential_symbol)
        return new(_check_atom_names(atom_names), _check_atom_separator(atom_separator),
                   _check_symbol(:identity, identity), _check_garside_symbol(garside_symbol),
                   _check_symbol(:bar_separator, bar_separator), _check_symbol(:empty_cell, empty_cell),
                   _check_symbol(:empty_cell_in_module, empty_cell_in_module),
                   _check_coefficient_style(coefficient_style),
                   _check_symbol(:differential_symbol, differential_symbol))
    end
end

_notation_error(msg) = throw(ArgumentError("Notation: " * msg))

# The checks of the settings, one method each with branches rather than one method per type:
# a method of one line that returns a constant would be folded away by the compiler, and
# the coverage would never see it.
function _check_atom_names(x)
    x === nothing && return nothing
    x isa Function && return x
    x === :subscripts && return x
    x isa AbstractVector ||
        _notation_error("atom_names = $(repr(x)), expected nothing, :subscripts, a vector of strings or a function")
    names = String[string(y) for y in x]
    for (a, name) in enumerate(names)
        isempty(name) && _notation_error("atom_names: the name of the atom $a is empty")
    end
    allunique(names) || _notation_error("atom_names: two atoms have the same name")
    return names
end

function _check_atom_separator(x)
    x === :auto && return x
    x isa AbstractString || _notation_error("atom_separator = $(repr(x)), expected :auto or a string")
    return String(x)
end

function _check_symbol(key::Symbol, x)
    x isa AbstractString || _notation_error("$key = $(repr(x)), expected a string")
    isempty(x) && _notation_error("$key is empty")
    return String(x)
end

_check_garside_symbol(x) = x === nothing ? nothing : _check_symbol(:garside_symbol, x)

const _COEFFICIENT_STYLES = (:auto, :times, :star, :parentheses)

function _check_coefficient_style(x)
    x isa Symbol && x in _COEFFICIENT_STYLES && return x
    _notation_error("coefficient_style = $(repr(x)), expected one of " * join((":$t" for t in _COEFFICIENT_STYLES), ", "))
end

# The defaults, in the order of the fields.
const _DEFAULT_NOTATION = Notation(nothing, :auto, "e", nothing, "|", "e", "[]", :auto, "d")

function Notation(N::Notation; settings...)
    for k in keys(settings)
        k in fieldnames(Notation) ||
            _notation_error("unknown setting $k; the settings are " * join(fieldnames(Notation), ", "))
    end
    return Notation((haskey(settings, k) ? settings[k] : getfield(N, k) for k in fieldnames(Notation))...)
end

Notation(; settings...) = Notation(_DEFAULT_NOTATION; settings...)

Base.:(==)(M::Notation, N::Notation) = all(isequal(getfield(M, k), getfield(N, k)) for k in fieldnames(Notation))
Base.hash(N::Notation, h::UInt) = foldr((k, a) -> hash(getfield(N, k), a), fieldnames(Notation); init = hash(:Notation, h))

# Only the settings that differ from the defaults.
function Base.show(io::IO, N::Notation)
    parts = String[]
    for k in fieldnames(Notation)
        v = getfield(N, k)
        isequal(v, getfield(_DEFAULT_NOTATION, k)) || push!(parts, "$k = $(repr(v))")
    end
    print(io, "Notation(", join(parts, ", "), ")")
end

# The settings of the session.
const _SESSION_NOTATION = Ref(_DEFAULT_NOTATION)

"""
    notation() -> Notation
    notation(io::IO) -> Notation

The settings of the display of the session ([`Notation`](@ref)), or those of the display
on `io`: the `Notation` under its key `:notation` if it has one, those of the session
otherwise.
"""
notation() = _SESSION_NOTATION[]

function notation(io::IO)
    N = get(io, :notation, _SESSION_NOTATION[])
    N isa Notation || throw(ArgumentError("the key :notation of the IOContext holds a $(typeof(N)), not a Notation"))
    return N
end

"""
    set_notation!(; settings...) -> Notation

Change the settings of the display of the session given as keywords (see
[`Notation`](@ref)), the others being kept; without argument, go back to the defaults.
Returns the new settings: `set_notation!(atom_names = ["a", "b"])` writes the atoms of a
germ with two atoms as letters, `set_notation!()` as numbers again.
"""
function set_notation!(; settings...)
    _SESSION_NOTATION[] = isempty(settings) ? _DEFAULT_NOTATION : Notation(_SESSION_NOTATION[]; settings...)
    return _SESSION_NOTATION[]
end

# ------------------------------------------------------------------------------------
# Names of the atoms, separators and markers, read by the displays and the readers of the
# monoids (for germs, in the layer of the resolution).

const _SUBSCRIPT_DIGITS = ('₀', '₁', '₂', '₃', '₄', '₅', '₆', '₇', '₈', '₉')

_subscript_string(n::Integer) = join(_SUBSCRIPT_DIGITS[d + 1] for d in reverse(digits(n)))

_atom_name(::Nothing, a::Int) = string(a)
_atom_name(::Symbol, a::Int) = "a" * _subscript_string(a)      # :subscripts, checked
_atom_name(names::Vector{String}, a::Int) = a <= length(names) ? names[a] : string(a)
function _atom_name(f::Function, a::Int)
    name = string(f(a))
    isempty(name) && _notation_error("atom_names: the function gives an empty name to the atom $a")
    return name
end

# The names of the atoms 1:n.
function _atom_names(N::Notation, n::Integer)
    names = String[_atom_name(N.atom_names, a) for a in 1:Int(n)]
    allunique(names) || _notation_error("atom_names: two atoms have the same name")
    return names
end

# Characters that may follow the first one in a name written side by side with others.
_is_modifier(c::Char) = c in _SUBSCRIPT_DIGITS || c in ('′', '″', '‴', '\'')

# Whether a name is one character, possibly followed by subscript digits or primes: names of
# this kind can be written side by side and split again before each first character.
function _is_self_delimiting(name::AbstractString)
    isempty(name) && return false
    _is_modifier(first(name)) && return false
    return all(_is_modifier, Iterators.drop(name, 1))
end

# The separator of the atoms of a word, for a germ whose atoms have the given names.
_atom_separator(N::Notation, names::Vector{String}) =
    N.atom_separator isa String ? N.atom_separator : (all(_is_self_delimiting, names) ? "" : "·")
