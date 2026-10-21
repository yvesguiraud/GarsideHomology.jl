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
