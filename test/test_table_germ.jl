# Tests of the tabulated germs (table.jl): agreement with the germs of explicit_germ.jl,
# built independently by brute force; independence of the provisional numbering; display;
# the checks of the constructor.

isdefined(Main, :ExplicitGerm) || include("explicit_germ.jl")

using Gauss
using Test

# Every answer of the interface of T is that of E, types included, for all simples.
function same_answers(E, T)
    S, N = Gauss.simples(E), natoms(E)
    Gauss.simples(T) == S && natoms(T) == N && nsimples(T) == nsimples(E) || return false
    for s in S
        length(T, s) === length(E, s) || return false
        Gauss.left_descents(T, s) == Gauss.left_descents(E, s) || return false
        Gauss.word(T, s) == Gauss.word(E, s) || return false
        all(Gauss.lcm_atom(T, s, a) === Gauss.lcm_atom(E, s, a) for a in 1:N) || return false
        for t in S
            Gauss.product(T, s, t) === Gauss.product(E, s, t) || return false
            Gauss.ldiv(T, s, t) === Gauss.ldiv(E, s, t) || return false
            Gauss.is_left_divisor(T, s, t) === Gauss.is_left_divisor(E, s, t) || return false
            Gauss.isless_lex(T, s, t) === Gauss.isless_lex(E, s, t) || return false
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
        T, renumber = Gauss.tabulate(E)
        @test T isa Gauss.TableGerm{max(1, cld(natoms(E), 64))}
        @test renumber == 0:(nsimples(E) - 1)
        @test same_answers(E, T)
        @test check_germ(T) == 0
        @test opposite(T) == Gauss.TableGerm(opposite(E))
        @test opposite(opposite(T)) == T
        @test Gauss.TableGerm(T) == T && hash(Gauss.TableGerm(T)) == hash(T)
        # Independence of the provisional numbering.
        R = right_covers_of(E)
        @test Gauss.TableGerm(R, 1) == T
        for seed in 1:3
            perm = _test_permutation(nsimples(E), seed)
            R2, identity2 = renumbered_covers(R, 1, perm)
            T2, renumber2 = Gauss.tabulate(R2, identity2)
            @test T2 == T
            @test all(renumber2[perm[p]] == renumber[p] for p in 1:nsimples(E))
        end
    end

    @testset "Sizes of the interval germs" begin
        B3 = Gauss.TableGerm(signed_interval_germ(3, :B, (-1, -2, -3)))
        @test nsimples(B3) == 20 && natoms(B3) == 9
        D4 = Gauss.TableGerm(signed_interval_germ(4, :D, (2, -4, -3, 1)))
        @test nsimples(D4) == 50 && natoms(D4) == 12         # Cat(D₄) = 50
        @test length(D4, nsimples(D4) - 1) == 4
    end

    @testset "Types and display" begin
        T = Gauss.TableGerm(classical_dihedral_germ(3))
        @test T isa Gauss.AbstractGerm && T isa Gauss.Complexes.AbstractMonoid
        @test sprint(show, T) == "TableGerm with 2 atoms and 6 simples"
        text = sprint(show, MIME("text/plain"), T)
        @test occursin("TableGerm{1} with 2 atoms and 6 simples", text)
        @test occursin("maximal length 3", text) && occursin("MiB", text)
        @test Gauss.TableGerm(dual_dihedral_germ(70)) isa Gauss.TableGerm{2}
        @test T != Gauss.TableGerm(dual_dihedral_germ(3))
        @test Gauss.product(T, 1, 4) === Int32(5) && Gauss.ldiv(T, 3, 5) === Int32(1)
        @test Gauss.product(T, 3, 4) === nothing && Gauss.ldiv(T, 4, 3) === nothing
        @test_throws BoundsError Gauss.product(T, 6, 1)
        @test_throws BoundsError Gauss.lcm_atom(T, 1, 3)
        A = MonoidAlgebra(Gauss.Complexes.AbstractAlgebra.ZZ, T)
        @test A.monoid === T && Gauss.Complexes.label_type(A) == NTuple{2, Int32}
    end

    @testset "Checks of the constructor" begin
        # Valid: the free monoid on one atom, truncated (simples e, a).
        G = Gauss.TableGerm([2 0], 1)
        @test natoms(G) == 1 && nsimples(G) == 2 && Gauss.lcm_atom(G, 0, 1) === (lcm = Int32(1), complement = Int32(1))
        @test_throws "no simple" Gauss.tabulate(zeros(Int32, 0, 0), 1)
        @test_throws "out of 0:2" Gauss.tabulate([3 0], 1)
        @test_throws "identity 3 out of 1:2" Gauss.tabulate([2 0], 3)
        @test_throws "atom 1 missing" Gauss.tabulate([0 2], 1)
        @test_throws "atom 1 is the identity" Gauss.tabulate([1 0], 1)
        @test_throws "atom 2 repeated" Gauss.tabulate([2 0; 2 0], 1)
        @test_throws "not reached" Gauss.tabulate([2 0 0], 1)
        @test_throws "not homogeneous" Gauss.tabulate([2 2], 1)
        # ab = bb.
        @test_throws "right cancellation" Gauss.tabulate([2 0 0 0; 3 4 4 0], 1)
        # aa = ab.
        @test_throws "left cancellation" Gauss.tabulate([2 4 0 0; 3 4 0 0], 1)
        # ab = ba = x and xa simple, aa not: a(ba) = (ab)a is simple, (aa)b is not.
        @test_throws "not associative" Gauss.tabulate([2 0 5 0 4; 3 5 0 0 0], 1)
        @test_throws "not bounded" Gauss.tabulate([2 0 0; 3 0 0], 1)
        chain(n) = reshape(Int32[(p < n ? p + 1 : 0) for p in 1:n], 1, n)
        @test length(Gauss.TableGerm(chain(256), 1), 255) == 255
        @test_throws "larger than 255" Gauss.tabulate(chain(257), 1)
        # The interval of a quasi-Coxeter element of D₄ is not a lattice: twelve pairs of
        # atoms have two minimal common multiples of length 3.
        S, atoms, e, mul, len = signed_interval(4, :D, (4, 3, -2, -1))
        @test length(S) == 54 && length(atoms) == 12
        R, identity = interval_covers(S, atoms, e, mul, len)
        @test_throws "no lcm" Gauss.tabulate(R, identity)
    end
end
