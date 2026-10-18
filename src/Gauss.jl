"""
    Gauss

Small free resolutions of Garside monoids and their homology.

Gauss computes the reduced lcm resolution of a Garside monoid, a variant of the
Dehornoy–Lafont resolution with far fewer cells, and the homology of the
corresponding groups with various coefficients.

This version is a package skeleton: the layers (germs, cell schemes, cell tree,
projection, resolution and certificates, modules, complexes) are added step by step.
"""
module Gauss

# Generic layer of free modules and chain complexes. It must not depend on anything
# else in Gauss, so that it can later become a package of its own.
include("complexes/Complexes.jl")

end # module Gauss
