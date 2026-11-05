# Working with other packages

```@meta
CurrentModule = GarsideHomology
```

## Typing the names

- `ℤ`: `\bbZ`, then Tab, in the Julia REPL and in the editors that know Julia (VS Code, Jupyter, Pluto, the Julia modes of Emacs and Vim).
- `ℤε`: `\bbZ`, Tab, `\varepsilon`, Tab. The other epsilon, `\epsilon`, gives `ℤϵ`, another character but the same object.
- `⊗`: `\otimes`, then Tab.

## Names shared with other packages

Some names of GarsideHomology are also exported by other packages, for other objects. They are public but not exported by GarsideHomology, or are to be qualified when these packages are loaded.

- `ZZ`, the ring of integers of AbstractAlgebra, is `ℤ` under its usual name. Nemo, Hecke and Oscar export a `ZZ` of their own: write `using GarsideHomology: ZZ` to use ours.
- `homology` and `differential` are exported by GarsideHomology, and also by Hecke and Oscar: with one of them loaded, write `GarsideHomology.homology` and `GarsideHomology.differential`, or import them explicitly (`using GarsideHomology: homology`).
- `tensor_product` (Hecke, Oscar) and `⊗` (Hecke, Oscar, Chevie): write `GarsideHomology.tensor_product`, or `using GarsideHomology: ⊗`.
- `basis`, `word`, `rank`, `elements` and `⊗`, exported by Chevie: with Chevie loaded, they are written without qualification, since the extension adds methods for our objects to the functions of Chevie. Without Chevie, write `using GarsideHomology: basis, word, rank, elements`. `elements` is also exported by FiniteFields, Hecke and Oscar.

## Limits

- Chevie represents the reflection groups by permutations of their roots, with 16-bit integers by default, and the extension of GarsideHomology stores the images of the roots with 16-bit integers too: the groups with more than 32 767 roots, such as G(29, 29, 7) and its 35 322 roots, are out of reach for now.
- A tabulated germ ([`TableGerm`](@ref)) takes about 16 bytes per simple and per atom (49 MB for B*(E₈)); above 5·10⁵ simples, `germ` gives a lazy germ for a classical braid monoid, and `dual_germ` an `ArgumentError`.
- A lazy germ ([`CoxeterGerm`](@ref)) numbers its simples as it meets them: it must not be queried from several threads at once.
- The products of the monoid beyond the simples, by normal forms, are not supported yet: `s * t` is defined only when it is a simple.
