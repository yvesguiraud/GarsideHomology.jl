# Display settings

```@meta
CurrentModule = GarsideHomology
```

The elements of the germs, the cells and the elements of the modules are written following settings that can be changed: the names of the atoms, the separators, the symbols of the identity, of Δ and of the differential, and the way the coefficients are written ([`Notation`](@ref)). The reading of written forms, `G("12")` or `P("2|12")`, follows the same settings. The examples continue those of the page [Resolutions, cells and elements](resolutions.md), in B⁺(A₂).

## Settings of the session

[`set_notation!`](@ref) changes the settings of the session given as keywords, and `set_notation!()` goes back to the defaults; [`notation`](@ref) reads them.

```julia
2 * d(P("2"))                    # 2×2[] - 2×[]2
set_notation!(atom_names = ["a", "b"])
2 * d(P("2"))                    # 2b[] - 2[]b
P("2b[] - 2[]b") == 2 * d(P("2"))   # true: the reading follows the settings
set_notation!()
```

With the atoms written as numbers, the marker `×` separates a coefficient from what follows; with letters, it is not needed.

## Settings of one display

A `Notation` under the key `:notation` of an `IOContext` applies to one display, and wins over the session:

```julia
show(IOContext(stdout, :notation => Notation(coefficient_style = :parentheses)), 2 * P(P("2")))   # 2([2])
show(IOContext(stdout, :notation => Notation(empty_cell_in_module = "[e]")), d(P("1")))           # 1[e] - [e]1
```

## The settings

- `atom_names`: the names of the atoms, by default their numbers; `:subscripts` for `a₁, a₂, …`, a vector of strings, or a function of the number of an atom.
- `atom_separator`: written between the atoms of a word. By default, the atoms are written side by side when each name is a single character, possibly followed by subscripts or primes, and separated by `·` otherwise, for instance in a germ with at least 10 atoms numbered: `1·12·3`.
- `identity`, `garside_symbol`: the identity (`e`), and Δ, written as a word by default or by a symbol such as `"Δ"`.
- `bar_separator`, `empty_cell`, `empty_cell_in_module`: the separator of the letters of a cell (`|`), the cell of degree 0 (`e`) and its generator in a module (`[]`).
- `coefficient_style`: `×` only when the names of the atoms can begin with a digit (the default), always `×` (`:times`) or `*` (`:star`), or parentheses around the monomial (`:parentheses`).
- `differential_symbol`: the name of the differential in the displays of the complexes, `d` by default (`"∂"`, for instance).
