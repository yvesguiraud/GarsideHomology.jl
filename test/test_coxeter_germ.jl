# Tests of the lazy germs of the classical braid monoids (coxeter.jl), without Chevie:
# the dihedral germs given by their Cartan matrices compared with ExplicitGerm, the
# numbering at the encounter, E₈, the errors.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using Gauss
using Test

const CARTAN_A2 = [2 -1; -1 2]
const CARTAN_B2 = [2 -2; -1 2]
const CARTAN_G2 = [2 -1; -3 2]
const CARTAN_E8 = [ 2  0 -1  0  0  0  0  0
                    0  2  0 -1  0  0  0  0
                   -1  0  2 -1  0  0  0  0
                    0 -1 -1  2 -1  0  0  0
                    0  0  0 -1  2 -1  0  0
                    0  0  0  0 -1  2 -1  0
                    0  0  0  0  0 -1  2 -1
                    0  0  0  0  0  0 -1  2]

@testset "CoxeterGerm" begin

    @testset "Dihedral germs by their Cartan matrices: $name" for (name, C, m) in
            [("A₂", CARTAN_A2, 3), ("B₂", CARTAN_B2, 4), ("B₂ transposed", permutedims(CARTAN_B2), 4),
             ("G₂", CARTAN_G2, 6)]
        G = Gauss.CoxeterGerm(C)
        @test G isa Gauss.CoxeterGerm{2} && G isa Gauss.AbstractGerm
        @test natoms(G) == 2 && nsimples(G) == 2m && Gauss.nnumbered(G) == 3
        @test Gauss.word(G, 1) == [1] && Gauss.word(G, 2) == [2] && length(G, 0) == 0
        @test same_germ(G, classical_dihedral_germ(m))
        @test Gauss.nnumbered(G) == 2m
        @test same_germ(opposite(G), classical_dihedral_germ(m))
    end

    @testset "E₈, lazily" begin
        G = Gauss.CoxeterGerm(CARTAN_E8)
        @test natoms(G) == 8 && nsimples(G) == 696729600 && Gauss.nnumbered(G) == 9
        Δ = foldl((ω, a) -> Gauss.lcm_atom(G, ω, a).lcm, 1:8; init = Int32(0))
        @test length(G, Δ) == 120 && collect(Gauss.left_descents(G, Δ)) == 1:8
        @test all(Gauss.product(G, Δ, a) === nothing for a in 1:8)
        @test Gauss.ldiv(G, Δ, Δ) === Int32(0) && Gauss.ldiv(G, 0, Δ) === Δ
        @test all(Gauss.lcm_atom(G, Δ, a) === (lcm = Δ, complement = Int32(0)) for a in 1:8)
        s = Gauss.product(G, Gauss.product(G, 1, 3)::Int32, 4)::Int32       # s₁s₃s₄
        @test length(G, s) == 3 && Gauss.ldiv(G, 1, s) == Gauss.product(G, 3, 4)
        @test Gauss.product(G, 3, 3) === nothing && Gauss.product(G, 1, 2) == Gauss.product(G, 2, 1)
        # A query numbers only the simple it returns.
        H = Gauss.CoxeterGerm(CARTAN_E8)
        walk(w) = foldl((x, a) -> Gauss.product(H, x, a)::Int32, w; init = Int32(0))
        s, t = walk([1, 3, 4]), walk([1, 3, 4, 5, 6])
        before = Gauss.nnumbered(H)
        @test Gauss.word(H, t) == [1, 3, 4, 5, 6] && Gauss.isless_lex(H, s, t) && !Gauss.isless_lex(H, t, s)
        @test Gauss.nnumbered(H) == before
        u = Gauss.ldiv(H, s, t)
        @test Gauss.nnumbered(H) == before + 1 && Gauss.word(H, u) == [5, 6]
        @test Gauss.ldiv(H, t, s) === nothing && Gauss.ldiv(H, walk([2, 4]), t) === nothing
        @test Gauss.nnumbered(H) == before + 2                  # s₂s₄ numbered by walk
        @test startswith(sprint(show, G), "CoxeterGerm of rank 8, ")
        @test G == G && G != Gauss.CoxeterGerm(CARTAN_E8)
        A = MonoidAlgebra(Gauss.Complexes.AbstractAlgebra.ZZ, G)
        @test A.monoid === G
    end

    @testset "Errors" begin
        @test_throws ArgumentError Gauss.CoxeterGerm([2 -1 0; -1 2 -1])           # not square
        @test_throws ArgumentError Gauss.CoxeterGerm([2 -1; -1 3])                 # diagonal
        @test_throws ArgumentError Gauss.CoxeterGerm([2 -1; 0 2])                  # not symmetric in its zeros
        @test_throws "not of finite crystallographic type" Gauss.CoxeterGerm([2 -2; -2 2])   # affine Ã₁
        @test_throws BoundsError Gauss.lcm_atom(Gauss.CoxeterGerm(CARTAN_A2), 0, 3)
    end
end
