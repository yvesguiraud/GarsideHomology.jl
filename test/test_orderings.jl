# Tests of the orderings (cells/orderings.jl): LCMOrdering and its order of the root, the
# preorders on the lcms, compatibility with divisibility, isless_deglex of the germs (also
# checked for every pair of simples by check_germ, explicit_germ.jl).

isdefined(Main, :check_germ) || include("explicit_germ.jl")

using GarsideHomology
using Test
import Chevie

@testset "Orderings" begin
    lex = FunctionPreorder((G, ω, x, y) -> GarsideHomology.isless_lex(G, x.complement, y.complement))

    @testset "LCMOrdering" begin
        S = LCMOrdering()
        @test S.preorder === ComplementDeglex() && S.root_order === nothing
        @test S == LCMOrdering(ComplementDeglex()) && hash(S) == hash(LCMOrdering(ComplementDeglex()))
        @test sprint(show, S) == "LCMOrdering(ComplementDeglex())"
        D = LCMOrdering(TrivialPreorder(); root_order = [2, 3, 1])
        @test D.root_order == Int32[2, 3, 1] && D.root_order isa Vector{Int32}
        @test sprint(show, D) == "LCMOrdering(TrivialPreorder(); root_order = [2, 3, 1])"
        @test D == LCMOrdering(TrivialPreorder(); root_order = Int32[2, 3, 1])
        @test hash(D) == hash(LCMOrdering(TrivialPreorder(); root_order = (2, 3, 1)))
        @test D != LCMOrdering(TrivialPreorder()) && D != LCMOrdering(TrivialPreorder(); root_order = 1:3)
        @test LCMOrdering(lex) == LCMOrdering(lex) && LCMOrdering(lex) != S
        p = [2, 3, 1]
        E = LCMOrdering(TrivialPreorder(); root_order = p)
        p[1] = 7
        @test E.root_order == Int32[2, 3, 1]                  # the order is copied
        @test_throws ArgumentError LCMOrdering(; root_order = [1, 1, 3])
        @test_throws ArgumentError LCMOrdering(; root_order = [0, 1])
        @test_throws ArgumentError LCMOrdering(; root_order = [1.5, 2])
        @test_throws MethodError LCMOrdering(:lex)
        G = dual_dihedral_germ(3)
        @test GarsideHomology._root_order(S, G) == Int32[1, 2, 3]
        @test GarsideHomology._root_order(D, G) == Int32[2, 3, 1]
        @test_throws ArgumentError GarsideHomology._root_order(LCMOrdering(; root_order = [2, 1]), G)
    end

    @testset "Preorders on B⁺(A₃)" begin
        G = germ(Chevie.BraidMonoid(Chevie.coxgroup(:A, 3)))
        x = GarsideHomology.lcm_atom(G, 1, 2)         # s₁ ∨ s₂ = s₁s₂s₁, complement s₂s₁
        y = GarsideHomology.lcm_atom(G, 1, 3)         # s₁ ∨ s₃ = s₁s₃, complement s₃
        @test GarsideHomology.word(G, x.complement) == [2, 1] && GarsideHomology.word(G, y.complement) == [3]
        @test GarsideHomology.precedes(ComplementDeglex(), G, 1, y, x)
        @test !GarsideHomology.precedes(ComplementDeglex(), G, 1, x, y)
        @test !GarsideHomology.precedes(ComplementDeglex(), G, 1, x, x)
        @test GarsideHomology.precedes(lex, G, 1, x, y) && !GarsideHomology.precedes(lex, G, 1, y, x)
        @test !GarsideHomology.precedes(TrivialPreorder(), G, 1, x, y)
        @test !GarsideHomology.precedes(TrivialPreorder(), G, 1, y, x)
        # No two lcms ω ∨ a, ω ∨ b of B⁺(A₃) are strictly comparable: every preorder is
        # compatible with divisibility, and the four variants of J0 coincide.
        @test GarsideHomology.extends_divisibility(ComplementDeglex(), G)
        @test GarsideHomology.extends_divisibility(TrivialPreorder(), G)
        @test GarsideHomology.extends_divisibility(lex, G)
        bad = FunctionPreorder((G, ω, x, y) -> 1)
        @test_throws TypeError GarsideHomology.precedes(bad, G, 1, x, y)
    end

    @testset "Preorders on B*(A₃)" begin
        G = GarsideHomology.TableGerm(Chevie.DualBraidMonoid(Chevie.coxgroup(:A, 3); c = collect(1:3)))
        @test GarsideHomology.extends_divisibility(ComplementDeglex(), G)
        @test !GarsideHomology.extends_divisibility(TrivialPreorder(), G)
        # The lexicographic order of the complements is not compatible: at ω = a₅,
        # a₅ ∨ a₁ ≺ a₅ ∨ a₄ = Δ, but the complement 2.6 of a₄ comes before the complement 6
        # of a₁. The construction never meets this pair (a₂ removes a₄ first), so that
        # Gauss_c and Gauss_a coincide on B*(A₃). Only ω = a₄ and ω = a₅ break the condition.
        @test !GarsideHomology.extends_divisibility(lex, G)
        x, y = GarsideHomology.lcm_atom(G, 5, 1), GarsideHomology.lcm_atom(G, 5, 4)
        @test GarsideHomology.word(G, x.complement) == [6] && GarsideHomology.word(G, y.complement) == [2, 6]
        @test x.lcm != y.lcm && GarsideHomology.is_left_divisor(G, x.lcm, y.lcm)
        @test GarsideHomology.precedes(lex, G, 5, y, x) && GarsideHomology.precedes(ComplementDeglex(), G, 5, x, y)
        @test !GarsideHomology.extends_divisibility(lex, G, [5]) && !GarsideHomology.extends_divisibility(lex, G, [4])
        @test GarsideHomology.extends_divisibility(lex, G, 0:3)
    end

    @testset "Preorders on explicit germs" begin
        D4 = signed_interval_germ(4, :D, (2, -4, -3, 1))      # B*(D₄)
        for G in (classical_dihedral_germ(4), dual_dihedral_germ(3), dual_dihedral_germ(5), D4)
            @test GarsideHomology.extends_divisibility(ComplementDeglex(), G)
        end
        # In B*(I₂(m)), the lcms of an atom with the others are all Δ: nothing to order.
        @test GarsideHomology.extends_divisibility(TrivialPreorder(), dual_dihedral_germ(5))
        @test !GarsideHomology.extends_divisibility(TrivialPreorder(), D4)
    end

    @testset "isless_deglex" begin
        T = GarsideHomology.TableGerm(Chevie.DualBraidMonoid(Chevie.coxgroup(:B, 3); c = collect(1:3)))
        generic(s, t) = invoke(GarsideHomology.isless_deglex, Tuple{GarsideHomology.AbstractGerm, Integer, Integer}, T, s, t)
        S = GarsideHomology.simples(T)
        @test all(GarsideHomology.isless_deglex(T, s, t) == generic(s, t) == (s < t) for s in S, t in S)
        @test_throws BoundsError GarsideHomology.isless_deglex(T, 0, nsimples(T))
        H = GarsideHomology.CoxeterGerm(Chevie.cartan(Chevie.coxgroup(:A, 3)))
        s12 = GarsideHomology.product(H, 1, 2)::Int32
        s21 = GarsideHomology.product(H, 2, 1)::Int32
        @test GarsideHomology.isless_deglex(H, 3, s12) && !GarsideHomology.isless_deglex(H, s12, 3)
        @test GarsideHomology.isless_deglex(H, s12, s21) && !GarsideHomology.isless_deglex(H, s21, s12)
        @test !GarsideHomology.isless_deglex(H, s12, s12)
    end

    @testset "Exports" begin
        for s in (:LCMOrdering, :ComplementDeglex, :TrivialPreorder, :FunctionPreorder)
            @test Base.isexported(GarsideHomology, s)
        end
        for s in (:AbstractLCMPreorder, :precedes, :extends_divisibility)
            @test Base.ispublic(GarsideHomology, s) && !Base.isexported(GarsideHomology, s)
        end
    end
end
