# Tests of the germs built from Chevie (extension GaussChevieExt): small classical and
# dual braid monoids compared completely with Chevie (chevie_germs.jl), consistency of
# the interface (check_germ of explicit_germ.jl), opposite monoids, the threshold of germ.

isdefined(Main, :compare_with_chevie) || include("chevie_germs.jl")
isdefined(Main, :check_germ) || include("explicit_germ.jl")

using Gauss
using Test
import Chevie

@testset "Germs from Chevie" begin
    A(n) = Chevie.coxgroup(:A, n)
    B(n) = Chevie.coxgroup(:B, n)
    cases = [
        "B⁺(A₃)" => Chevie.BraidMonoid(A(3)),
        "B⁺(B₃)" => Chevie.BraidMonoid(B(3)),
        "B⁺(H₃)" => Chevie.BraidMonoid(Chevie.coxgroup(:H, 3)),
        "B*(A₃)" => Chevie.DualBraidMonoid(A(3); c = collect(1:3)),
        "B*(B₃)" => Chevie.DualBraidMonoid(B(3); c = collect(1:3)),
        "B*(G₄)" => Chevie.DualBraidMonoid(Chevie.crg(4)),
        "B*(G₅)" => Chevie.DualBraidMonoid(Chevie.crg(5)),
        "B*(G₆)" => Chevie.DualBraidMonoid(Chevie.crg(6)),
    ]
    sizes = Dict("B⁺(A₃)" => 24, "B⁺(B₃)" => 48, "B⁺(H₃)" => 120, "B*(A₃)" => 14,
                 "B*(B₃)" => 20, "B*(G₄)" => 5, "B*(G₅)" => 6, "B*(G₆)" => 8)

    @testset "$name" for (name, M) in cases
        G = Gauss.TableGerm(M)
        @test G isa Gauss.TableGerm
        @test nsimples(G) == sizes[name]
        @test germ(M) == G
        r = compare_with_chevie(G, M)
        @test r.failures == 0 && r.complete
        @test check_germ(G) == 0
    end

    @testset "Enumeration without Chevie's simples (c3bis): $name" for (name, M) in [
            "B*(A₃)" => Chevie.DualBraidMonoid(A(3); c = collect(1:3)),
            "B*(B₃)" => Chevie.DualBraidMonoid(B(3); c = collect(1:3)),
            "B*(D₄)" => Chevie.DualBraidMonoid(Chevie.coxgroup(:D, 4); c = collect(1:4)),
            "B*(H₃)" => Chevie.DualBraidMonoid(Chevie.coxgroup(:H, 3); c = collect(1:3)),
            "B*(I₂(5))" => Chevie.DualBraidMonoid(Chevie.coxgroup(:I, 2, 5); c = [1, 2]),
            "B*(G₄)" => Chevie.DualBraidMonoid(Chevie.crg(4)),
            "B*(G₅)" => Chevie.DualBraidMonoid(Chevie.crg(5)),
            "B*(G₆)" => Chevie.DualBraidMonoid(Chevie.crg(6)),        # c: a word of 19 letters
            "B*(G₂₅)" => Chevie.DualBraidMonoid(Chevie.crg(25)),
            "B*(G₂₆)" => Chevie.DualBraidMonoid(Chevie.crg(26))]
        ext = Base.get_extension(Gauss, :GaussChevieExt)
        G = Gauss.TableGerm(M)
        # Same germ as by the simples of left_divisors (the reference of c3).
        @test G == Gauss.TableGerm(ext._right_covers(M; method = :chevie)...)
        # From the group alone, for the same Coxeter element.
        @test dual_germ(M.W; c = M.δword) == G
    end

    @testset "Default Coxeter element of dual_germ" begin
        # Chevie's default c is 1:n for these groups (j1/c3bis_coxeter.jl): Springer's test
        # accepts it and the germ is that of DualBraidMonoid(W).
        for W in (A(3), B(3), Chevie.coxgroup(:H, 3), Chevie.crg(4), Chevie.crg(5), Chevie.crg(25))
            @test dual_germ(W) == Gauss.TableGerm(Chevie.DualBraidMonoid(W; c = collect(1:length(Chevie.gens(W)))))
        end
        @test_throws "well-generated" dual_germ(Chevie.crg(31))     # 5 generators in rank 4
        # Cat(A₁₁) = 208 012 simples, decided from the degrees.
        @test_throws "above the tabulation threshold" dual_germ(Chevie.coxgroup(:A, 11))
    end

    @testset "Opposite monoids" begin
        # w ↦ w⁻¹ fixes the atoms: B⁺(W)ᵒᵖ ≅ B⁺(W) and, for a real group,
        # B*(W, c)ᵒᵖ ≅ B*(W, c⁻¹), with the same order of the atoms.
        W = B(3)
        @test opposite(Gauss.TableGerm(Chevie.BraidMonoid(W))) == Gauss.TableGerm(Chevie.BraidMonoid(W))
        @test opposite(Gauss.TableGerm(Chevie.DualBraidMonoid(W; c = [1, 2, 3]))) ==
              Gauss.TableGerm(Chevie.DualBraidMonoid(W; c = [3, 2, 1]))
        @test opposite(Gauss.TableGerm(Chevie.DualBraidMonoid(W; c = [1, 2, 3]))) !=
              Gauss.TableGerm(Chevie.DualBraidMonoid(W; c = [1, 2, 3]))
    end

    @testset "Entry point" begin
        # Above the threshold, without enumeration: 696 729 600 simples.
        @test_throws "above the tabulation threshold" germ(Chevie.BraidMonoid(Chevie.coxgroup(:E, 8)))
        @test_throws MethodError germ(A(2))
        ext = Base.get_extension(Gauss, :GaussChevieExt)
        # (B*(G₃₄) is not built here: Chevie takes minutes; j1/c3_chevie.jl of CatLab
        # compares the numbers of simples with Chevie on every case.)
        @test ext._simples_count(Chevie.DualBraidMonoid(Chevie.crg(37))) == 25080
        @test ext._simples_count(Chevie.BraidMonoid(Chevie.coxgroup(:E, 8))) == 696729600
        # Timings of the phases, for the measurements of c3.
        timings = Pair{Symbol, Float64}[]
        R, identity = ext._right_covers(Chevie.DualBraidMonoid(A(3); c = collect(1:3)), timings)
        G, _ = Gauss._tabulate(R, identity, timings)
        @test first.(timings) == [:enumeration, :covers, :input, :lengths, :right_cancellation,
                                  :left_covers, :left_quotients, :bounded, :canonical,
                                  :lex_ranks, :tables, :lcms]
        @test all(>=(0), last.(timings))
        @test G == Gauss.TableGerm(Chevie.DualBraidMonoid(A(3); c = collect(1:3)))
    end
end
