# Tests of the elements of a germ for the user (src/germs/elements.jl, J1.h): construction by
# calling the germ, operations, simples and elements, writing with the settings of the
# display, reading, the Garside element; the sign module ℤε; the methods of the functions of
# Chevie of the same names (package extension).

using GarsideHomology
using GarsideHomology: GermElement, word, elements, CoxeterGerm
using Test
import Chevie

# A germ that numbers none of its simples and does not count them (for the checks).
struct _BareGerm <: GarsideHomology.AbstractGerm end

@testset "Elements of germs" begin
    A2 = germ(Chevie.coxgroup(:A, 2))                    # simples e, 1, 2, 12, 21, 121
    C2 = CoxeterGerm([2 -1; -1 2])                       # the same monoid, lazy
    D65 = GarsideHomology._dual_dihedral_germ(65)        # B*(I₂(65)): 65 atoms, c = 1·2

    @testset "Construction" begin
        @test isone(A2()) && A2() == one(A2(1)) && A2(1).index == 1 && A2(1).germ === A2
        @test A2(1, 2).index == 3 && A2(2, 1).index == 4 && A2(1, 2, 1) == A2(2, 1, 2) == GermElement(A2, 5)
        @test_throws ArgumentError A2(1, 1)                   # not a simple
        @test_throws ArgumentError A2(1, 3)                   # no atom 3
        @test_throws ArgumentError A2(0)
        @test_throws ArgumentError A2(7)
        @test occursin("G(1, 2)", sprint(showerror, try A2(12) catch e; e end))
        @test occursin("G(\"121\")", sprint(showerror, try A2(121) catch e; e end))
        @test occursin("G(\"1·2·1\")", sprint(showerror, try D65(121) catch e; e end))
        @test_throws ArgumentError GermElement(A2, 6)
        @test_throws ArgumentError GermElement(A2, -1)
        @test GermElement(C2, 0) == C2()
        @test_throws ArgumentError GermElement(C2, 100)       # beyond the simples numbered so far
        @test GermElement(_BareGerm(), 5).index == 5          # nothing to check against
        @test_throws ArgumentError GermElement(_BareGerm(), -1)
    end

    @testset "Operations" begin
        @test word(A2(1, 2, 1)) == Int32[1, 2, 1] && length(A2(1, 2, 1)) == 3 && word(A2()) == Int32[]
        @test A2(1) * A2(2) == A2(1, 2) && A2() * A2(2) == A2(2) && A2(2) * A2() == A2(2)
        @test_throws ArgumentError A2(1) * A2(1)
        @test A2(1) \ A2(1, 2, 1) == A2(2, 1)
        @test_throws ArgumentError A2(2) \ A2(1, 2)
        @test_throws ArgumentError A2(1) * C2(1)              # different germs
        @test A2(1) != C2(1) && hash(A2(1, 2)) == hash(A2(1, 2)) && length(Set([A2(1), A2(1), A2(2)])) == 2
        @test sort([A2(1, 2, 1), A2(2, 1), A2(2), A2(), A2(1, 2), A2(1)]) == collect(simples(A2))
        @test length.(simples(A2)) == [0, 1, 1, 2, 2, 3]
        @test [A2(1)] .* A2(2) == [A2(1, 2)]                  # a scalar for broadcasting
        @test C2(1) * C2(2) * C2(1) == C2(1, 2, 1) && C2(1) \ C2(1, 2) == C2(2)
    end

    @testset "Simples and elements" begin
        S = simples(A2)
        @test S isa AbstractVector{<:GermElement} && length(S) == 6 && S[end] == A2(1, 2, 1)
        @test string.(S) == ["e", "1", "2", "12", "21", "121"]
        @test elements(A2) == S
        @test startswith(sprint(summary, S), "6-element simples of TableGerm")
        @test_throws ArgumentError simples(C2)                # a lazy germ does not list its simples
        @test_throws BoundsError S[7]
    end

    @testset "Writing" begin
        show_with(N, x) = sprint(show, x; context = :notation => N)
        @test string(A2(1, 2, 1)) == "121" && string(A2()) == "e" && repr(A2(2, 1)) == "21"
        @test string(D65(1, 2)) == "1·2" && string(D65(12)) == "12"
        @test show_with(Notation(atom_names = ["a", "b"]), A2(1, 2, 1)) == "aba"
        @test show_with(Notation(atom_names = :subscripts), A2(1, 2)) == "a₁a₂"
        @test show_with(Notation(atom_separator = "."), A2(1, 2, 1)) == "1.2.1"
        @test show_with(Notation(identity = "1"), A2()) == "1"
        Δ = Notation(garside_symbol = "Δ")
        @test show_with(Δ, A2(1, 2, 1)) == "Δ" && show_with(Δ, A2(1, 2)) == "12"
        @test show_with(Δ, C2(2, 1, 2)) == "Δ" && show_with(Δ, D65(1, 2)) == "Δ"
        @test GarsideHomology._garside_element(C2) == C2(1, 2, 1).index
        @test invoke(GarsideHomology._garside_element, Tuple{GarsideHomology.AbstractGerm}, D65) == 66
    end

    @testset "Reading" begin
        @test A2("121") == A2("1·2·1") == A2("1.2.1") == A2(" 1 2 1 ") == A2("(121)") == A2("1⋅2•1") == A2(1, 2, 1)
        @test A2("e") == A2() && A2("12") == A2(1, 2) && A2("((21))") == A2(2, 1)
        @test D65("1·2") == D65("1.2") == D65(1, 2) && D65("12") == D65(12)
        @test C2("212") == C2(1, 2, 1)
        for bad in ("11", "3", "", "()", "x", "(1)(2)", "((1)", "1·3")
            @test_throws ArgumentError A2(bad)
        end
        @test_throws ArgumentError D65("121")                 # the atom 121, which D65 has not
        try
            set_notation!(atom_names = ["a", "b"], garside_symbol = "Δ")
            @test A2("aba") == A2("a·b·a") == A2(1, 2, 1) && A2("ba") == A2(2, 1) && A2("Δ") == A2(1, 2, 1)
            @test A2("12") == A2(1, 2)                         # numbers are read too
            @test string(A2(1, 2)) == "ab"                     # the settings of the session
            set_notation!(atom_separator = "-")
            @test A2("a-b-a") == A2(1, 2, 1) && string(A2(1, 2)) == "a-b"
        finally
            set_notation!()
        end
        @test all(A2(string(s)) == s for s in simples(A2))
        @test all(D65(string(s)) == s for s in simples(D65))
    end

    @testset "Sign module ℤε" begin
        @test ℤε === ℤϵ && ℤε isa SignModule && GarsideHomology.base_ring(ℤε) === ℤ
        R = lcm_resolution(A2)
        @test homology(R, ℤε) == homology(R, SignModule())
        for s in (:ℤε, :ℤϵ)
            @test Base.isexported(GarsideHomology, s)
        end
    end

    @testset "Functions of Chevie of the same names" begin
        @test Chevie.word === Chevie.PermGroups.word && Chevie.elements === Chevie.PermGroups.elements
        @test Chevie.basis === Chevie.Algebras.basis && Chevie.rank === Chevie.Symbols.rank
        @test Chevie.word(A2(1, 2)) == Int32[1, 2] && Chevie.word(A2, 4) == Int32[2, 1]
        @test Chevie.elements(A2) == simples(A2)
        R = lcm_resolution(A2)
        @test Chevie.basis(R[1]) == GarsideHomology.basis(R[1])
        @test Chevie.rank(R[1]) == 2 && Chevie.rank(homology(R, ℤ)[1]) == 1 && Chevie.rank(ℤε) == 1
    end

    @testset "Exports" begin
        for s in (:GermElement, :elements)
            @test Base.ispublic(GarsideHomology, s) && !Base.isexported(GarsideHomology, s)
        end
    end
end
