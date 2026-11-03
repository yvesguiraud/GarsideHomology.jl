# Tests of the detailed displays (src/complexes/display.jl, J1.h, point 5b): complexes and
# resolutions with their cells and boundaries, differentials, long elements alone and in
# vectors; the complete view, and the view cut to the screen (:limit).

using GarsideHomology
using GarsideHomology: basis, BasedFreeModuleElement
using Test
import Chevie

const _DI = GarsideHomology.Complexes
_limited(rows, columns) = (:limit => true, :displaysize => (rows, columns))

@testset "Detailed displays" begin
    A2 = germ(Chevie.coxgroup(:A, 2))                    # B⁺(A₂): cells e; 1, 2; 2|12
    P = lcm_resolution(A2)
    d = differential(P)
    x = P("2|12")
    cells = "\nCells:\n  0: e\n  1: 1, 2\n  2: 2|12"
    top = "\n  d(2|12) = " * string(d(x))

    @testset "Resolutions" begin
        @test sprint(show, MIME"text/plain"(), P) ==
              sprint(show, P) * cells * "\nDifferential:\n  d(1) = 1[] - []1\n  d(2) = 2[] - []2" * top
        @test sprint(show, MIME"text/plain"(), P; context = _limited(24, 200)) ==
              sprint(show, P) * cells * "\nDifferential:\n  (boundaries of degree 1 omitted, such as d(1) = 1[] - []1)" * top
        @test occursin("\n  ∂(2|12) = ", sprint(show, MIME"text/plain"(), P; context = :notation => Notation(differential_symbol = "∂")))
        D65 = lcm_resolution(GarsideHomology._dual_dihedral_germ(65))     # cells 1, 65, 64
        lines = split(sprint(show, MIME"text/plain"(), D65; context = _limited(12, 80)), '\n')
        @test length(lines) == 10 && lines[end - 1] == "  ⋮ (62 more boundaries)"
        @test startswith(lines[7], "  (boundaries of degree 1 omitted, such as d(1) = ") && occursin(", …, ", lines[4])
        @test startswith(lines[8], "  d(") && startswith(lines[10], "  d(")
        narrow = split(sprint(show, MIME"text/plain"(), D65; context = _limited(40, 30)), '\n')
        @test all(textwidth(l) <= 30 for l in narrow if startswith(l, "  d("))
        @test length(split(sprint(show, MIME"text/plain"(), D65), '\n')) == 1 + 1 + 3 + 1 + 65 + 64
    end

    @testset "Differentials" begin
        @test sprint(show, MIME"text/plain"(), d) ==
              "Differential of " * sprint(show, P) * "\n  d(1) = 1[] - []1\n  d(2) = 2[] - []2" * top
        @test sprint(show, MIME"text/plain"(), d; context = _limited(24, 200)) ==
              "Differential of " * sprint(show, P) * "\n  (boundaries of degree 1 omitted, such as d(1) = 1[] - []1)" * top
        E = extend_scalars(P)
        @test occursin("such as d(1) = 0)", sprint(show, MIME"text/plain"(), E; context = _limited(24, 200)))
    end

    @testset "Long elements" begin
        D65 = lcm_resolution(GarsideHomology._dual_dihedral_germ(65))
        v = sum(D65(c) for c in basis(D65[1]))               # [1] + [2] + ⋯ + [65]
        @test v isa BasedFreeModuleElement && length(v) == 65
        @test sprint(show, v; context = _limited(24, 30)) == "[1] + [2] + [3] … (65 terms)"
        plain = sprint(show, MIME"text/plain"(), v; context = _limited(6, 30))
        @test startswith(plain, "[1] + [2]") && endswith(plain, "[8] … (65 terms)")
        @test count('+', sprint(show, MIME"text/plain"(), v)) == 64 && count('+', string(v)) == 64
        @test count("… (65 terms)", sprint(show, MIME"text/plain"(), [v, v]; context = _limited(24, 30))) == 2
        @test sprint(show, zero(v); context = _limited(24, 30)) == "0"
        @test sprint(show, D65(basis(D65[1])[1]); context = _limited(24, 2)) == "[1]"   # the first term, always
    end
end
