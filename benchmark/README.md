# Benchmarks

Measurement harnesses and protocols for GarsideHomology.jl; they are not part of the loaded package.

Protocol: one fresh Julia session per measurement, at least three runs in random order,
minimum kept; machine, Julia version and commit written in every CSV line; results saved
incrementally. Timings are compared only on the same machine and Julia version.

## Construction of the lcm resolution (J1.d, task d7)

`construction.jl` measures one round of one case in a fresh session, `run_construction.jl`
runs the campaign: the seven reference cases of J0 (Gauss_a), two paths each (P1 through
the monoid of Chevie, as J0; P2 with the package alone), three rounds in random order,
then the minima beside those of J0. Results in `benchmark/out/` (not in git). The columns
of the CSV are in English since J1.g (`case`, `path`, `round`, `cells`, `terms`, `chains`,
`commit`, `dirty`, then `t_<step>_s`, `gc_<step>_s`, `compile_<step>_s` and the memory in
GiB for the steps `group`, `monoid`, `germ`, `tree`, `diff`): a CSV written before, with
French columns (d7 on N9), is not resumed. The CSV of J0 keeps its own columns.

```
julia --project=benchmark -e "using Pkg; Pkg.instantiate()"
julia --project=benchmark benchmark/run_construction.jl N9
```
