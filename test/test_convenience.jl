# Tests of the conveniences of J1.g (task g1): germ(W) of a finite Coxeter group of Chevie,
# lcm_resolution and homology of a monoid of Chevie through its germ, homology of a germ,
# ncells of a resolution.

isdefined(Main, :check_germ) || include("explicit_germ.jl")

using GarsideHomology
using GarsideHomology: CoxeterGerm, cell_tree, ncells
using Test
import Chevie

@testset "Conveniences" begin
    G = FinitelyGeneratedAbelianGroup
    GG = FinitelyGeneratedGradedAbelianGroup
    trivial3 = GG([G(1), G(1), G(0)])                 # H_*(B₃; ℤ), B₃ the braid group on three strands
    sign3 = GG([G(0, [2]), G(0, [3]), G(0)])          # H_*(B₃; ℤε)

    @testset "germ of a finite Coxeter group" begin
        for W in (Chevie.coxgroup(:A, 3), Chevie.coxgroup(:B, 3))
            @test germ(W) == germ(Chevie.BraidMonoid(W))
        end
        @test germ(Chevie.coxgroup(:E, 8)) isa CoxeterGerm          # 696 729 600 simples: lazy
        @test germ(Chevie.crg(23)) == germ(Chevie.coxgroup(:H, 3))  # crg(23) is coxgroup(:H, 3) in Chevie
        @test_throws "see dual_germ(W)" germ(Chevie.crg(4))
        @test_throws ArgumentError germ(Chevie.crg(3, 3, 3))       # a PermRootGroup
    end

    @testset "lcm_resolution of a monoid of Chevie" begin
        for M in (Chevie.BraidMonoid(Chevie.coxgroup(:A, 3)), Chevie.DualBraidMonoid(Chevie.coxgroup(:B, 3)),
                  Chevie.DualBraidMonoid(Chevie.crg(4)))
            for ordering in (LCMOrdering(), LCMOrdering(TrivialPreorder()))
                R, S = lcm_resolution(M, ordering), lcm_resolution(germ(M), ordering)
                @test cell_tree(R).ordering == ordering
                @test range(R) == range(S)
                @test all(ncells(R, n) == ncells(cell_tree(S), n) for n in range(R))
                @test all(differential(R, n) == differential(S, n) for n in first(range(R)):(last(range(R)) - 1))
            end
        end
        @test lcm_resolution(Chevie.BraidMonoid(Chevie.coxgroup(:A, 2))) isa GarsideHomology.LCMResolution
    end

    @testset "ncells of a resolution" begin
        R = lcm_resolution(dual_germ(Chevie.crg(4)))                 # B*(G₄): cells 1, 3, 2
        @test ncells(R) == 6 == ncells(cell_tree(R)) && [ncells(R, n) for n in range(R)] == [1, 3, 2]
        @test ncells(R, 3) == 0 && ncells(R, -1) == 0                 # outside the degrees of R
        D = lcm_resolution(Chevie.BraidMonoid(Chevie.coxgroup(:A, 3)), LCMOrdering(TrivialPreorder()))
        @test ncells(D) == sum(ncells(D, n) for n in range(D)) && ncells(D, 0) == 1
    end

    @testset "homology of a germ and of a monoid" begin
        # B₃ as B⁺(A₂), B*(A₂) and B*(G₄), by germs and by monoids of Chevie.
        for X in (dual_dihedral_germ(3), classical_dihedral_germ(3), dual_germ(Chevie.crg(4)),
                  Chevie.BraidMonoid(Chevie.coxgroup(:A, 2)), Chevie.DualBraidMonoid(Chevie.crg(4)))
            @test homology(X) == trivial3 == homology(X, TrivialModule())
            @test homology(X, SignModule()) == sign3
            @test homology(X, SignModule(); value_type = Int64) == sign3     # keywords to the solver
            @test homology(X, 1) == G(1) && homology(X, SignModule(), 1) == G(0, [3])
            @test homology(X; ordering = LCMOrdering(TrivialPreorder())) == trivial3
            @test homology(X, SignModule(), 0; ordering = LCMOrdering(TrivialPreorder())) == G(0, [2])
        end
        # H_*(B₄; ℤ) by the classical braid monoid of A₃, tabulated and lazy.
        trivial4 = GG([G(1), G(1), G(0, [2]), G(0)])
        @test homology(germ(Chevie.coxgroup(:A, 3))) == trivial4 == homology(CoxeterGerm(Chevie.cartan(Chevie.coxgroup(:A, 3))))
    end
end
