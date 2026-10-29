# Tests of the tabulated germs (table.jl): agreement with the germs of explicit_germ.jl,
# built independently by brute force; independence of the provisional numbering; display;
# the checks of the constructor.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using GarsideHomology
using GarsideHomology: MonoidAlgebra    # public, not exported (Oscar exports a MonoidAlgebra)
using Test

# Every answer of the interface of T is that of E, types included, for all simples.
function same_answers(E, T)
    S, N = GarsideHomology.simples(E), natoms(E)
    GarsideHomology.simples(T) == S && natoms(T) == N && nsimples(T) == nsimples(E) || return false
    for s in S
        length(T, s) === length(E, s) || return false
        GarsideHomology.left_descents(T, s) == GarsideHomology.left_descents(E, s) || return false
        GarsideHomology.word(T, s) == GarsideHomology.word(E, s) || return false
        all(GarsideHomology.lcm_atom(T, s, a) === GarsideHomology.lcm_atom(E, s, a) for a in 1:N) || return false
        for t in S
            GarsideHomology.product(T, s, t) === GarsideHomology.product(E, s, t) || return false
            GarsideHomology.ldiv(T, s, t) === GarsideHomology.ldiv(E, s, t) || return false
            GarsideHomology.is_left_divisor(T, s, t) === GarsideHomology.is_left_divisor(E, s, t) || return false
            GarsideHomology.isless_lex(T, s, t) === GarsideHomology.isless_lex(E, s, t) || return false
        end
    end
    return true
end

# The covers R renumbered by the permutation `perm` of 1:n (p becomes perm[p]).
function renumbered_covers(R, identity, perm)
    N, n = size(R)
    R2 = zeros(Int32, N, n)
    for p in 1:n, a in 1:N
        q = R[a, p]
        R2[a, perm[p]] = q == 0 ? 0 : perm[q]
    end
    return R2, perm[identity]
end

# A pseudo-random permutation of 1:n, reproducible from its seed.
function _test_permutation(n, seed)
    x = UInt64(seed)
    values = [(x = x * 0x5851f42d4c957f2d + 0x14057b7ef767814f) for _ in 1:n]
    return sortperm(values)
end

const TABLE_TEST_GERMS = [
    "B⁺(A₂)" => classical_dihedral_germ(3), "B*(A₂)" => dual_dihedral_germ(3),
    "B⁺(I₂(4))" => classical_dihedral_germ(4), "B*(I₂(4))" => dual_dihedral_germ(4),
    "B⁺(I₂(5))" => classical_dihedral_germ(5), "B*(I₂(5))" => dual_dihedral_germ(5),
    "B⁺(I₂(70))" => classical_dihedral_germ(70), "B*(I₂(70))" => dual_dihedral_germ(70),
    "[1, -1] in B₃" => signed_interval_germ(3, :B, (-1, -2, -3)),
    "B*(D₄)" => signed_interval_germ(4, :D, (2, -4, -3, 1)),
]

@testset "TableGerm" begin

    @testset "Agreement with ExplicitGerm: $name" for (name, E) in TABLE_TEST_GERMS
        T, renumber = GarsideHomology.tabulate(E)
        @test T isa GarsideHomology.TableGerm{max(1, cld(natoms(E), 64))}
        @test renumber == 0:(nsimples(E) - 1)
        @test same_answers(E, T)
        @test check_germ(T) == 0
        @test opposite(T) == GarsideHomology.TableGerm(opposite(E))
        @test opposite(opposite(T)) == T
        @test GarsideHomology.TableGerm(T) == T && hash(GarsideHomology.TableGerm(T)) == hash(T)
        # Independence of the provisional numbering.
        R = right_covers_of(E)
        @test GarsideHomology.TableGerm(R, 1) == T
        for seed in 1:3
            perm = _test_permutation(nsimples(E), seed)
            R2, identity2 = renumbered_covers(R, 1, perm)
            T2, renumber2 = GarsideHomology.tabulate(R2, identity2)
            @test T2 == T
            @test all(renumber2[perm[p]] == renumber[p] for p in 1:nsimples(E))
        end
    end

    @testset "Sizes of the interval germs" begin
        B3 = GarsideHomology.TableGerm(signed_interval_germ(3, :B, (-1, -2, -3)))
        @test nsimples(B3) == 20 && natoms(B3) == 9
        D4 = GarsideHomology.TableGerm(signed_interval_germ(4, :D, (2, -4, -3, 1)))
        @test nsimples(D4) == 50 && natoms(D4) == 12         # Cat(D₄) = 50
        @test length(D4, nsimples(D4) - 1) == 4
    end

    @testset "Types and display" begin
        T = GarsideHomology.TableGerm(classical_dihedral_germ(3))
        @test T isa GarsideHomology.AbstractGerm && T isa GarsideHomology.Complexes.AbstractMonoid
        @test sprint(show, T) == "TableGerm with 2 atoms and 6 simples"
        text = sprint(show, MIME("text/plain"), T)
        @test occursin("TableGerm{1} with 2 atoms and 6 simples", text)
        @test occursin("maximal length 3", text) && occursin("MiB", text)
        @test GarsideHomology.TableGerm(dual_dihedral_germ(70)) isa GarsideHomology.TableGerm{2}
        @test T != GarsideHomology.TableGerm(dual_dihedral_germ(3))
        @test GarsideHomology.product(T, 1, 4) === Int32(5) && GarsideHomology.ldiv(T, 3, 5) === Int32(1)
        @test GarsideHomology.product(T, 3, 4) === nothing && GarsideHomology.ldiv(T, 4, 3) === nothing
        @test_throws BoundsError GarsideHomology.product(T, 6, 1)
        @test_throws BoundsError GarsideHomology.lcm_atom(T, 1, 3)
        A = MonoidAlgebra(GarsideHomology.Complexes.AbstractAlgebra.ZZ, T)
        @test A.monoid === T && GarsideHomology.Complexes.label_type(A) == NTuple{2, Int32}
    end

    @testset "Checks of the constructor" begin
        # Valid: the free monoid on one atom, truncated (simples e, a).
        G = GarsideHomology.TableGerm([2 0], 1)
        @test natoms(G) == 1 && nsimples(G) == 2 && GarsideHomology.lcm_atom(G, 0, 1) === (lcm = Int32(1), complement = Int32(1))
        @test_throws "no simple" GarsideHomology.tabulate(zeros(Int32, 0, 0), 1)
        @test_throws "out of 0:2" GarsideHomology.tabulate([3 0], 1)
        @test_throws "identity 3 out of 1:2" GarsideHomology.tabulate([2 0], 3)
        @test_throws "atom 1 missing" GarsideHomology.tabulate([0 2], 1)
        @test_throws "atom 1 is the identity" GarsideHomology.tabulate([1 0], 1)
        @test_throws "atom 2 repeated" GarsideHomology.tabulate([2 0; 2 0], 1)
        @test_throws "not reached" GarsideHomology.tabulate([2 0 0], 1)
        @test_throws "not homogeneous" GarsideHomology.tabulate([2 2], 1)
        # ab = bb.
        @test_throws "right cancellation" GarsideHomology.tabulate([2 0 0 0; 3 4 4 0], 1)
        # aa = ab.
        @test_throws "left cancellation" GarsideHomology.tabulate([2 4 0 0; 3 4 0 0], 1)
        # ab = ba = x and xa simple, aa not: a(ba) = (ab)a is simple, (aa)b is not.
        @test_throws "not associative" GarsideHomology.tabulate([2 0 5 0 4; 3 5 0 0 0], 1)
        @test_throws "not bounded" GarsideHomology.tabulate([2 0 0; 3 0 0], 1)
        chain(n) = reshape(Int32[(p < n ? p + 1 : 0) for p in 1:n], 1, n)
        @test length(GarsideHomology.TableGerm(chain(256), 1), 255) == 255
        @test_throws "larger than 255" GarsideHomology.tabulate(chain(257), 1)
        # The interval of a quasi-Coxeter element of D₄ is not a lattice: twelve pairs of
        # atoms have two minimal common multiples of length 3.
        S, atoms, e, mul, len = signed_interval(4, :D, (4, 3, -2, -1))
        @test length(S) == 54 && length(atoms) == 12
        R, identity = interval_covers(S, atoms, e, mul, len)
        @test_throws "no lcm" GarsideHomology.tabulate(R, identity)
    end
end
