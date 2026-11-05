# Chevie extension

```@meta
CurrentModule = GarsideHomologyChevieExt
```

Loading [Chevie.jl](https://github.com/jmichel7/Chevie.jl), by `using Chevie`, loads this extension. The methods it adds are documented with the functions they extend: [`germ`](@ref GarsideHomology.germ), [`TableGerm`](@ref GarsideHomology.TableGerm), [`dual_germ`](@ref GarsideHomology.dual_germ), [`lcm_resolution`](@ref GarsideHomology.lcm_resolution) and [`homology`](@ref GarsideHomology.homology) for the groups and monoids of Chevie, and the methods of the functions `⊗`, `word`, `elements`, `basis` and `rank` of Chevie for the objects of GarsideHomology. The docstring of the extension describes how the simples are enumerated.

```@docs
GarsideHomologyChevieExt
```
