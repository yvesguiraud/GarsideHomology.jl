# Orderings, certificates and snapshots

```@meta
CurrentModule = GarsideHomology
```

## Orderings

An [`LCMOrdering`](@ref) fixes the cells of an lcm resolution: at each lcm met by the construction, a preorder on the lcms with the atoms decides which atoms give cells, ties keeping the order of the atoms at the root.

```julia
R = lcm_resolution(G)                                         # the reduced resolution
D = lcm_resolution(G, LCMOrdering(TrivialPreorder()))         # Dehornoy–Lafont, order 1:N
D′ = lcm_resolution(G, LCMOrdering(TrivialPreorder(); root_order = p))   # for the order p
```

- [`ComplementDeglex`](@ref), the default, orders the lcms by their complements, the shorter first, then by smallest word; it gives the reduced resolution.
- [`TrivialPreorder`](@ref) ties all the lcms: the construction then keeps the order of the root everywhere, and gives the order resolution of Dehornoy and Lafont for that order.
- [`FunctionPreorder`](@ref) takes any function, to explore other preorders. A preorder that is not compatible with divisibility may give a resolution that is not reduced, which is to be checked by its certificate; [`extends_divisibility`](@ref) tests the compatibility.

`ncells(R)` and `ncells(R, n)` give the numbers of cells, to compare two orderings.

## Certificates

[`certify`](@ref)`(R)` checks the resolution and returns a [`ResolutionCertificate`](@ref), shown as a short report:
- correctness: homogeneity of the boundaries, ∂² = 0 computed in the bimodules with the products of the germ, and Euler characteristic of the cells of each weight against the Möbius function μ(1, ω);
- minimality, the absence of terms `1·y·1` in the differential, and linearity, every term `σ y τ` having `ℓ(σ) + ℓ(τ) = 1`.

[`is_correct`](@ref), [`is_minimal`](@ref) and [`is_linear`](@ref) read the report. A failed check is reported, never thrown. The check of ∂² = 0 is the costliest: `certify(R; square_zero = false)` skips it, and `nworkers` sets the number of tasks it uses. The Möbius function is computed by the crosscut theorem of Rota up to 16 atoms, and by a pass over the intervals above (keyword `mobius`).

## Snapshots

[`write_snapshot`](@ref)`(path, R)` writes a resolution in a compressed text format (`path` ending in `.txt.gz`), with its cells and its differential, and by default its fingerprints beside it (`.sha.csv`). This is the format of the reference resolutions of the tests; it is described in the file `test/data/reference/FORMAT.md` of the [repository](https://github.com/yvesguiraud/GarsideHomology.jl).
