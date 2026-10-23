# Tests of the germ layer: bit sets (bitsets.jl), the interface of germs and its generic
# methods (interface.jl), on the test germ ExplicitGerm of explicit_germ.jl, built for the
# classical and dual braid monoids of I₂(m); exports.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using Gauss
using Gauss: MonoidAlgebra    # public, not exported (Oscar exports a MonoidAlgebra)
using Test
import Chevie

const SBS = Gauss.StaticBitSet

# A small deterministic generator (xorshift), so that a failure can be reproduced.
mutable struct _XorShift
    x::UInt64
end
function _next!(g::_XorShift)
    x = g.x
    x ⊻= x << 13
    x ⊻= x >> 7
    x ⊻= x << 17
    g.x = x
    return x
end
_random_subset(g, n, p) = [i for i in 1:n if _next!(g) % 1000 < p]

_sum_elements(s) = (c = 0; for i in s; c += i; end; c)

@testset "Germs" begin

    @testset "StaticBitSet" begin
        @test isbitstype(SBS{2}) && sizeof(SBS{2}) == 16
        @test SBS{2} <: AbstractSet{Int}
        @test_throws ArgumentError SBS{0}()
        @test_throws ArgumentError SBS{1}([65])
        @test_throws ArgumentError SBS{1}([0])
        @test_throws ArgumentError SBS{1}([1.5])
        @test_throws ArgumentError first(SBS{1}())
        @test isempty(SBS{3}()) && length(SBS{3}()) == 0 && collect(SBS{3}()) == Int[]
        s = SBS{2}([65, 1, 64, 128, 1])
        @test collect(s) == [1, 64, 65, 128] && length(s) == 4 && first(s) == 1
        @test 64 in s && 65 in s && !(63 in s) && !(0 in s) && !(129 in s) && !(-1 in s)
        @test sprint(show, s) == "StaticBitSet{2}([1, 64, 65, 128])"
        @test s == Set([1, 64, 65, 128]) && s != SBS{2}([1, 64, 65])
        @test copy(s) === s
        @test _sum_elements(s) == 258
        @test (@allocated _sum_elements(s)) == 0

        g = _XorShift(0x9e3779b97f4a7c15)
        for K in 1:3, trial in 1:40
            n = 64K
            p = (50, 300, 700, 1000)[mod1(trial, 4)]
            a, b = _random_subset(g, n, p), _random_subset(g, n, 1000 - p)
            trial == 1 && (a = [1, 64, n]; b = [n])                  # word boundaries
            A, B = BitSet(a), BitSet(b)
            sa, sb = SBS{K}(a), SBS{K}(b)
            @test collect(sa) == collect(A) && length(sa) == length(A) && isempty(sa) == isempty(A)
            isempty(A) || @test first(sa) == first(A)
            @test all((i in sa) == (i in A) for i in -1:(n + 2))
            @test collect(sa ∪ sb) == collect(A ∪ B)
            @test collect(sa ∩ sb) == collect(A ∩ B)
            @test collect(setdiff(sa, sb)) == collect(setdiff(A, B))
            @test collect(symdiff(sa, sb)) == collect(symdiff(A, B))
            @test (sa ⊆ sb) == (A ⊆ B) && (sa ∩ sb ⊆ sa) && isdisjoint(sa, sb) == isdisjoint(A, B)
            @test union(sa, sb, sa) == sa ∪ sb && intersect(sa, sb, sb) == sa ∩ sb
            @test (sa == sb) == (A == B) && sa == SBS{K}(collect(sa))
            @test hash(sa) == hash(SBS{K}(reverse(collect(sa))))
        end
    end

    @testset "B⁺(A₂), written by hand" begin
        G = classical_dihedral_germ(3)           # e, s, t, st, ts, sts = 0, …, 5
        @test G isa Gauss.AbstractGerm && G isa Gauss.Complexes.AbstractMonoid
        @test natoms(G) == 2 && nsimples(G) == 6 && simples(G) == 0:5
        @test Gauss.atoms(G) == 1:2 && Gauss.atoms(G) isa UnitRange{Int32}
        @test one(G) === Int32(0) && isone(G, 0) && !isone(G, 3)
        @test [Gauss.word(G, s) for s in 0:5] == [Int32[], [1], [2], [1, 2], [2, 1], [1, 2, 1]]
        @test Gauss.word(G, 5) isa Vector{Int32}
        @test [length(G, s) for s in 0:5] == [0, 1, 1, 2, 2, 3]
        @test sort(0:5; lt = (s, t) -> Gauss.isless_lex(G, s, t)) == [0, 1, 3, 5, 2, 4]
        @test Gauss.product(G, 1, 2) === Int32(3) && Gauss.product(G, 2, 1) === Int32(4)
        @test Gauss.product(G, 3, 1) === Int32(5) && Gauss.product(G, 1, 4) === Int32(5)
        @test Gauss.product(G, 2, 3) === Int32(5) && Gauss.product(G, 4, 2) === Int32(5)
        @test Gauss.product(G, 1, 1) === nothing && Gauss.product(G, 3, 2) === nothing
        @test Gauss.product(G, 0, 5) === Int32(5) && Gauss.product(G, 5, 0) === Int32(5)
        @test Gauss.ldiv(G, 1, 5) === Int32(4) && Gauss.ldiv(G, 2, 5) === Int32(3)
        @test Gauss.ldiv(G, 3, 5) === Int32(1) && Gauss.ldiv(G, 4, 5) === Int32(2)
        @test Gauss.ldiv(G, 1, 4) === nothing && Gauss.ldiv(G, 5, 5) === Int32(0)
        @test Gauss.is_left_divisor(G, 1, 3) && !Gauss.is_left_divisor(G, 2, 3)
        @test Gauss.is_left_divisor(G, 3, 5) && !Gauss.is_left_divisor(G, 3, 4)
        @test !Gauss.is_left_divisor(G, 5, 3) && Gauss.is_left_divisor(G, 0, 4)
        @test collect(Gauss.left_descents(G, 5)) == [1, 2] && collect(Gauss.left_descents(G, 3)) == [1]
        @test isempty(Gauss.left_descents(G, 0)) && Gauss.left_descents(G, 4) isa SBS{1}
        @test Gauss.lcm_atom(G, 0, 2) == (lcm = 2, complement = 2)
        @test Gauss.lcm_atom(G, 1, 2) == (lcm = 5, complement = 4)
        @test Gauss.lcm_atom(G, 3, 1) == (lcm = 3, complement = 0)
        @test Gauss.lcm_atom(G, 3, 2) == (lcm = 5, complement = 1)
        @test Gauss.lcm_atom(G, 4, 1) == (lcm = 5, complement = 2)
        ωa, τ = Gauss.lcm_atom(G, 1, 2)
        @test ωa === Int32(5) && τ === Int32(4)
        # w ↦ w⁻¹ is an anti-automorphism fixing the atoms: the opposite germ is the same.
        @test opposite(G) == G
        @test check_germ(G) == 0
        # A germ describes the monoid of a monoid algebra; its labels are Int32.
        A = MonoidAlgebra(Gauss.Complexes.AbstractAlgebra.ZZ, G, :right)
        @test Gauss.Complexes.label_type(A) == Int32 && A.monoid === G
    end

    @testset "B*(A₂), written by hand" begin
        G = dual_dihedral_germ(3)                # e, a, b, c, Δ = 0, …, 4; ab = bc = ca = Δ
        @test natoms(G) == 3 && nsimples(G) == 5
        @test [Gauss.word(G, s) for s in 0:4] == [Int32[], [1], [2], [3], [1, 2]]
        @test sort(0:4; lt = (s, t) -> Gauss.isless_lex(G, s, t)) == [0, 1, 4, 2, 3]
        @test Gauss.product(G, 1, 2) === Int32(4) && Gauss.product(G, 2, 3) === Int32(4)
        @test Gauss.product(G, 3, 1) === Int32(4)
        @test Gauss.product(G, 2, 1) === nothing && Gauss.product(G, 1, 3) === nothing
        @test Gauss.product(G, 1, 1) === nothing && Gauss.product(G, 4, 1) === nothing
        @test Gauss.ldiv(G, 1, 4) === Int32(2) && Gauss.ldiv(G, 2, 4) === Int32(3)
        @test Gauss.ldiv(G, 3, 4) === Int32(1) && Gauss.ldiv(G, 1, 2) === nothing
        @test collect(Gauss.left_descents(G, 4)) == [1, 2, 3]
        # The example of the discussion of I2: atom 3 is c, the complement 2 is b.
        @test Gauss.lcm_atom(G, 1, 3) == (lcm = 4, complement = 2)
        @test Gauss.lcm_atom(G, 2, 1) == (lcm = 4, complement = 3)
        @test Gauss.lcm_atom(G, 4, 2) == (lcm = 4, complement = 0)
        @test Gauss.lcm_atom(G, 0, 3) == (lcm = 3, complement = 3)
        # Opposite: ba = cb = ac = Δ, so that a ·ᵒᵖ c = ca = Δ and the word of Δ is 1.3.
        O = opposite(G)
        @test O != G && opposite(O) == G
        @test Gauss.product(O, 1, 3) === Int32(4) && Gauss.product(O, 1, 2) === nothing
        @test Gauss.word(O, 4) == [1, 3]
        @test check_germ(G) == 0 && check_germ(O) == 0
    end

    @testset "Dihedral germs, m = $m" for m in (4, 5, 70)
        B = classical_dihedral_germ(m)
        @test natoms(B) == 2 && nsimples(B) == 2m
        @test length(B, nsimples(B) - 1) == m && opposite(B) == B
        @test check_germ(B) == 0
        D = dual_dihedral_germ(m)
        @test natoms(D) == m && nsimples(D) == m + 2      # Cat(I₂(m)) = m + 2
        @test Gauss.left_descents(D, m + 1) isa SBS{cld(m, 64)}
        @test collect(Gauss.left_descents(D, m + 1)) == 1:m
        @test Gauss.lcm_atom(D, 1, m) == (lcm = m + 1, complement = 2)
        @test Gauss.product(D, m, 1) === Int32(m + 1)      # rₘr₁ = Δ
        @test check_germ(D) == 0
    end

    @testset "Hurwitz orbit" begin
        # The transpositions 1 = (12), 2 = (13), 3 = (23) of S₃; xyx⁻¹ for x ≠ y is the
        # third one. The reduced decompositions of a 3-cycle: 3 = n!·hⁿ/|W| = 2·9/6.
        conj(i, j) = i == j ? i : 6 - i - j
        orbit = Gauss.hurwitz_orbit((1, 3), conj)
        @test orbit[1] == (1, 3) && Set(orbit) == Set([(1, 3), (2, 1), (3, 2)])
        @test Gauss.hurwitz_orbit((Int32(2),), conj) == [(Int32(2),)]
        @test_throws ArgumentError Gauss.hurwitz_orbit((1, 3), conj; limit = 2)
    end

    @testset "Exports" begin
        for s in (:natoms, :nsimples, :simples, :opposite, :germ, :dual_germ)
            @test Base.isexported(Gauss, s)
        end
        for s in (:AbstractGerm, :StaticBitSet, :product, :ldiv, :left_descents, :lcm_atom,
                  :atoms, :is_left_divisor, :word, :isless_lex, :isless_deglex, :TableGerm, :tabulate,
                  :hurwitz_orbit, :CoxeterGerm, :nnumbered, :verify, :count_multichains, :GermFingerprint, :fingerprint)
            @test Base.ispublic(Gauss, s) && !Base.isexported(Gauss, s)
        end
        # No name exported by Gauss is exported by Chevie for another object.
        collisions = Symbol[]
        for s in names(Gauss)
            (s === :Gauss || !Base.isexported(Gauss, s)) && continue
            Base.isexported(Chevie, s) && isdefined(Chevie, s) || continue
            getglobal(Chevie, s) === getglobal(Gauss, s) || push!(collisions, s)
        end
        @test isempty(collisions)
        # Names exported by Oscar, Hecke or Chevie for other objects, which Gauss keeps public
        # without exporting them (decisions of J1.b to J1.e; the sources of those packages
        # are checked by hand for each new name). Exported on purpose despite a collision:
        # homology (Oscar) and differential (Hecke, and Oscar, which reexports Hecke).
        for s in (:MonoidAlgebra, :tensor, :tensor_product, :⊗, :verify, :product, :atoms, :word,
                  :rank, :trivial_module)
            @test !Base.isexported(Gauss, s)
        end
    end
end
