# Tests of the lazy germs of the classical braid monoids (coxeter.jl), without Chevie:
# the dihedral germs given by their Cartan matrices compared with ExplicitGerm, the
# numbering at the encounter, E₈, the errors.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using GarsideHomology
using GarsideHomology: MonoidAlgebra    # public, not exported (Oscar exports a MonoidAlgebra)
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
        G = GarsideHomology.CoxeterGerm(C)
        @test G isa GarsideHomology.CoxeterGerm{2} && G isa GarsideHomology.AbstractGerm
        @test natoms(G) == 2 && nsimples(G) == 2m && GarsideHomology.nnumbered(G) == 3
        @test GarsideHomology.word(G, 1) == [1] && GarsideHomology.word(G, 2) == [2] && length(G, 0) == 0
        @test same_germ(G, classical_dihedral_germ(m))
        @test GarsideHomology.nnumbered(G) == 2m
        @test same_germ(opposite(G), classical_dihedral_germ(m))
    end

    @testset "E₈, lazily" begin
        G = GarsideHomology.CoxeterGerm(CARTAN_E8)
        @test natoms(G) == 8 && nsimples(G) == 696729600 && GarsideHomology.nnumbered(G) == 9
        Δ = foldl((ω, a) -> GarsideHomology.lcm_atom(G, ω, a).lcm, 1:8; init = Int32(0))
        @test length(G, Δ) == 120 && collect(GarsideHomology.left_descents(G, Δ)) == 1:8
        @test all(GarsideHomology.product(G, Δ, a) === nothing for a in 1:8)
        @test GarsideHomology.ldiv(G, Δ, Δ) === Int32(0) && GarsideHomology.ldiv(G, 0, Δ) === Δ
        @test all(GarsideHomology.lcm_atom(G, Δ, a) === (lcm = Δ, complement = Int32(0)) for a in 1:8)
        s = GarsideHomology.product(G, GarsideHomology.product(G, 1, 3)::Int32, 4)::Int32       # s₁s₃s₄
        @test length(G, s) == 3 && GarsideHomology.ldiv(G, 1, s) == GarsideHomology.product(G, 3, 4)
        @test GarsideHomology.product(G, 3, 3) === nothing && GarsideHomology.product(G, 1, 2) == GarsideHomology.product(G, 2, 1)
        # A query numbers only the simple it returns.
        H = GarsideHomology.CoxeterGerm(CARTAN_E8)
        walk(w) = foldl((x, a) -> GarsideHomology.product(H, x, a)::Int32, w; init = Int32(0))
        s, t = walk([1, 3, 4]), walk([1, 3, 4, 5, 6])
        before = GarsideHomology.nnumbered(H)
        @test GarsideHomology.word(H, t) == [1, 3, 4, 5, 6] && GarsideHomology.isless_lex(H, s, t) && !GarsideHomology.isless_lex(H, t, s)
        @test GarsideHomology.nnumbered(H) == before
        u = GarsideHomology.ldiv(H, s, t)
        @test GarsideHomology.nnumbered(H) == before + 1 && GarsideHomology.word(H, u) == [5, 6]
        @test GarsideHomology.ldiv(H, t, s) === nothing && GarsideHomology.ldiv(H, walk([2, 4]), t) === nothing
        @test GarsideHomology.nnumbered(H) == before + 2                  # s₂s₄ numbered by walk
        @test startswith(sprint(show, G), "CoxeterGerm of rank 8, ")
        @test G == G && G != GarsideHomology.CoxeterGerm(CARTAN_E8)
        A = MonoidAlgebra(GarsideHomology.Complexes.AbstractAlgebra.ZZ, G)
        @test A.monoid === G
    end

    @testset "Errors" begin
        @test_throws ArgumentError GarsideHomology.CoxeterGerm([2 -1 0; -1 2 -1])           # not square
        @test_throws ArgumentError GarsideHomology.CoxeterGerm([2 -1; -1 3])                 # diagonal
        @test_throws ArgumentError GarsideHomology.CoxeterGerm([2 -1; 0 2])                  # not symmetric in its zeros
        @test_throws "not of finite crystallographic type" GarsideHomology.CoxeterGerm([2 -2; -2 2])   # affine Ã₁
        @test_throws BoundsError GarsideHomology.lcm_atom(GarsideHomology.CoxeterGerm(CARTAN_A2), 0, 3)
        A2 = GarsideHomology.CoxeterGerm(CARTAN_A2)
        @test GarsideHomology._first_descent(A2.keys[1], A2.npositive) == 0           # the identity
    end
end
