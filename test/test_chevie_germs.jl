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
        # For G₂₉, 1.2.3.4 has order h = 20 but is not regular: the search finds 1.2.4.3,
        # Chevie's choice.
        W29 = Chevie.crg(29)
        G29 = @test_logs (:info, r"not a Coxeter element.*\[1, 2, 4, 3\]") dual_germ(W29)
        @test G29 == Gauss.TableGerm(Chevie.DualBraidMonoid(W29; c = [1, 2, 4, 3]))
        ext = Base.get_extension(Gauss, :GaussChevieExt)
        p = [1, 3, 2]
        @test ext._next_permutation!(p) && p == [2, 1, 3]
        p = [3, 2, 1]
        @test !ext._next_permutation!(p)
        # No order of the generators of A₁ × A₁ is regular for e^{2iπ/2} = −1: the eigenspace
        # of s₁s₂ = −1 is the whole plane.
        @test_throws "no order of the generators" ext._default_coxeter_word(A(1) * A(1), Val(2))
        # −1 = c³ in G₄ has no fixed vector but is not a product of two reflections (r⁻¹·(−1)
        # never has a fixed vector): the search of a reduced decomposition fails.
        W4 = Chevie.crg(4)
        @test_throws "no decomposition" ext._reduced_decomposition(W4, ext._Roots(W4, 2), W4(1, 2)^3,
                                                                   Int[], unique(Chevie.refls(W4)), Val(2))
        @test_throws "well-generated" dual_germ(Chevie.crg(31))     # 5 generators in rank 4
        # Cat(A₁₂) = 742 900 simples, decided from the degrees: no lazy dual germ.
        @test_throws "above the tabulation threshold" dual_germ(Chevie.coxgroup(:A, 12))
    end

    @testset "CoxeterGerm against TableGerm: $name" for (name, W) in [
            "B⁺(A₃)" => A(3), "B⁺(B₃)" => B(3), "B⁺(D₄)" => Chevie.coxgroup(:D, 4),
            "B⁺(A₄)" => A(4), "B⁺(B₄)" => B(4), "B⁺(F₄)" => Chevie.coxgroup(:F, 4)]
        L = Gauss.CoxeterGerm(Chevie.cartan(W))
        @test same_germ(L, Gauss.TableGerm(Chevie.BraidMonoid(W)))
    end

    @testset "verify and Fuss–Catalan numbers with Chevie" begin
        ext = Base.get_extension(Gauss, :GaussChevieExt)
        for W in (A(3), B(3), Chevie.coxgroup(:H, 3), Chevie.crg(4), Chevie.crg(25))
            G = dual_germ(W; c = W isa Chevie.FiniteCoxeterGroup ? collect(1:length(Chevie.gens(W))) : nothing)
            @test Gauss.verify(G; full = true) === nothing
            @test all(Gauss.count_multichains(G, k) == ext._fuss_catalan(W, k) for k in 1:3)
        end
        @test [Gauss.count_multichains(dual_germ(A(3)), k) for k in 1:3] == [14, 55, 140]
        @test Gauss.verify(Gauss.TableGerm(Chevie.BraidMonoid(B(3))); full = true) === nothing
    end

    @testset "Fingerprints with Chevie" begin
        fp = Gauss.fingerprint
        # B*(Bₙ) by coxgroup(:B, n) and by crg(2, 1, n): two presentations, one germ.
        for n in 2:5
            @test fp(dual_germ(B(n))) == fp(dual_germ(Chevie.crg(2, 1, n)))
        end
        @test fp(dual_germ(B(3); c = [1, 2, 3])) == fp(dual_germ(B(3); c = [3, 2, 1]))
        fA3, fB3 = fp(dual_germ(A(3))), fp(dual_germ(B(3)))
        @test fA3 != fp(Gauss.TableGerm(Chevie.BraidMonoid(A(3)))) && fA3 != fB3
        # B*(B₃): μ(1, Δ) = (−1)³ Π (dᵢ + h − 2)/dᵢ = −10; [1, −1] in B₃ has −7.
        @test fB3.mobius == -10 && fB3.natoms == 9 && sum(fB3.simples_by_length) == 20
    end

    @testset "Isomorphisms of dual germs (f4)" begin
        iso = Gauss.germ_isomorphism
        B3 = dual_germ(B(3))
        # G(m, 1, n): the dual germ of Bₙ, through the path of Hurwitz (2 is not a degree).
        @test iso(dual_germ(Chevie.crg(3, 1, 3)), B3) !== nothing
        @test iso(dual_germ(Chevie.crg(4, 1, 3)), B3) !== nothing
        @test iso(dual_germ(Chevie.crg(3, 1, 2)), dual_germ(B(2))) !== nothing
        # crg(2, 1, n) is the root datum of SO(2n + 1), crg(2, 2, n) that of SO(2n).
        @test iso(dual_germ(Chevie.crg(2, 1, 3)), B3) !== nothing
        @test iso(dual_germ(Chevie.crg(2, 2, 4)), dual_germ(Chevie.coxgroup(:D, 4))) !== nothing
        @test iso(dual_germ(Chevie.crg(2, 2, 3)), dual_germ(A(3))) !== nothing
        # G(3, 3, 3): 18 simples and 8 atoms, against 20 and 9 for B₃.
        @test iso(dual_germ(Chevie.crg(3, 3, 3)), B3) === nothing
    end

    @testset "Root datum of GL_n: crg(1, 1, n)" begin
        # Degrees 1, 2, …, n: the Fuss–Catalan numbers are those of the type Aₙ₋₁.
        G = dual_germ(Chevie.crg(1, 1, 3))
        @test nsimples(G) == 5 && natoms(G) == 3
        @test Gauss.germ_isomorphism(G, dual_germ(A(2))) !== nothing
        @test nsimples(dual_germ(Chevie.crg(1, 1, 4))) == 14
        ext = Base.get_extension(Gauss, :GaussChevieExt)
        @test [ext._fuss_catalan(Chevie.crg(1, 1, 4), k) for k in 1:3] == [14, 55, 140]
    end

    @testset "Orders of the Weyl groups from the Cartan matrices" begin
        for (T, n) in [(:A, 12), (:B, 10), (:D, 10), (:E, 6), (:E, 7), (:E, 8), (:F, 4), (:G, 2)]
            W = Chevie.coxgroup(T, n)
            @test nsimples(Gauss.CoxeterGerm(Chevie.cartan(W))) == prod(Chevie.degrees(W))
        end
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
        # Above the threshold, a braid monoid gets a lazy CoxeterGerm, without enumeration.
        E8 = germ(Chevie.BraidMonoid(Chevie.coxgroup(:E, 8)))
        @test E8 isa Gauss.CoxeterGerm{8} && nsimples(E8) == 696729600 && Gauss.nnumbered(E8) == 9
        @test germ(Chevie.BraidMonoid(Chevie.coxgroup(:B, 7))) isa Gauss.CoxeterGerm    # 645 120 simples
        @test Gauss._TABULATION_THRESHOLD == 5 * 10^5
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
