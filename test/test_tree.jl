# Tests of the tree of the cells (cells/tree.jl) and of the C lines of the snapshots
# (snapshots.jl): structure of the tree, C lines written by hand for B⁺(A₂) and B*(A₂),
# and the fingerprints of the C lines of small cases of J0, for the four variants
# (strict comparison for Gauss_a, Gauss_b, Gauss_d; in bar form for Gauss_c).

isdefined(Main, :check_germ) || include("explicit_germ.jl")
isdefined(Main, :read_fingerprints) || include("snapshot_io.jl")

using Gauss
using Gauss: CellTree, ncells     # public, not exported
using Test
import Chevie

const TREE_REFERENCE = joinpath(@__DIR__, "data", "reference")

# The four variants of J0 as orderings: Gauss_a and Gauss_d by the public preorders,
# Gauss_b and Gauss_c by preorders of the tests (length, then lexicographic order of the
# complements), with the form of the comparison (decisions of 22 October 2026): strict for
# a and d, in bar form for c, numbers of cells by degree for b, whose variant of J0 combines
# two orders. For the classical monoids, Gauss_a is checked by the invariant of
# minimality instead (check_classical): its children are ordered by deglex, those of J0 by
# lex, which differ from E₆ on.
const TREE_VARIANTS = [
    "a" => (LCMOrdering(ComplementDeglex()), :atoms),
    "b" => (LCMOrdering(FunctionPreorder((G, ω, x, y) -> length(G, x.complement) < length(G, y.complement))), :counts),
    "c" => (LCMOrdering(FunctionPreorder((G, ω, x, y) -> Gauss.isless_lex(G, x.complement, y.complement))), :bar),
    "d" => (LCMOrdering(TrivialPreorder()), :atoms),
]

# The fingerprints of the C lines of T, degree => (lines, atoms, bar).
function tree_c_fingerprints(T)
    out = Dict{Int, Tuple{Int, String, String}}()
    for n in 0:Gauss._top_degree(T)
        a = sort!(Gauss._cell_lines(T, n))
        b = sort!(Gauss._cell_lines(T, n; bar = true))
        out[n] = (length(a), _sha_lines(a), _sha_lines(b))
    end
    return out
end

function reference_c_fingerprints(case, variant)
    rows, _ = read_fingerprints(joinpath(TREE_REFERENCE, "$(case)_$(variant).sha.csv"))
    return Dict(parse(Int, r.n) => (r.lines, r.atoms, r.bar) for r in rows if r.section == "C")
end

# The invariant of minimality of a classical braid monoid with N atoms: one cell for each
# subset J of the atoms (weight w₀(J), degree |J|), hence 2ᴺ cells, binomial(N, k) in
# degree k, and pairwise distinct weights.
function check_classical(T)
    N = natoms(T.germ)
    return ncells(T) == 2^N && all(ncells(T, k) == binomial(N, k) for k in 0:N) &&
           allunique(Gauss.cell_weight(T, x) for x in 1:ncells(T))
end

# The structure of a tree: contiguous children, parents, degrees, lists read in the tree,
# letters and weights given by the lcms, bar forms and atom sequences.
function check_tree(T)
    failures = 0
    fail(msg) = (failures += 1; @warn msg; nothing)
    G = T.germ
    n = ncells(T)
    sum(Gauss._degree_counts(T)) == n || fail("degree counts")
    vcat([collect(Gauss.cell_range(T, d)) for d in 0:Gauss._top_degree(T)]...) == 1:n || fail("cell ranges")
    Gauss.cell_degree(T, 1) == 0 && Gauss.cell_weight(T, 1) == 0 || fail("root")
    next = 2
    for x in 1:n
        ch = Gauss.cell_children(T, x)
        isempty(ch) || first(ch) == next || fail("children of $x not contiguous")
        next += length(ch)
        for (k, y) in enumerate(ch)
            T.parent[y] == x || fail("parent of $y")
            Gauss.cell_degree(T, y) == Gauss.cell_degree(T, x) + 1 || fail("degree of $y")
            collect(Gauss.admissible_atoms(T, y)) == [Gauss.cell_atom(T, z) for z in ch[1:(k - 1)]] ||
                fail("list of $y")
            Gauss.cell_atom(T, y) in Gauss.admissible_atoms(T, x) || fail("atom of $y not in the list of $x")
            r = Gauss.lcm_atom(G, Gauss.cell_weight(T, x), Gauss.cell_atom(T, y))
            (r.lcm, r.complement) == (Gauss.cell_weight(T, y), Gauss.cell_letter(T, y)) ||
                fail("letter and weight of $y")
        end
        b = Gauss.bar_form(T, x)
        length(b) == Gauss.cell_degree(T, x) == length(Gauss.atom_sequence(T, x)) || fail("bar form of $x")
        foldl((s, t) -> Gauss.product(G, s, t)::Int32, b; init = Int32(0)) == Gauss.cell_weight(T, x) ||
            fail("product of the letters of $x")
    end
    next == n + 1 || fail("children: $next")
    return failures
end

@testset "Cell tree" begin
    @testset "B⁺(A₂) and B*(A₂), written by hand" begin
        # The C lines of the tiny reference snapshots classic_A2_a and dual_A2_a.
        T = CellTree(classical_dihedral_germ(3))
        @test sort!(vcat([Gauss._cell_lines(T, n) for n in 0:2]...)) ==
              ["C 0 - [] e", "C 1 1 [1] 1", "C 1 2 [2] 2", "C 2 2.1 [2|1.2] 1.2.1"]
        D = CellTree(dual_dihedral_germ(3))
        @test sort!(vcat([Gauss._cell_lines(D, n) for n in 0:2]...)) ==
              ["C 0 - [] e", "C 1 1 [1] 1", "C 1 2 [2] 2", "C 1 3 [3] 3", "C 2 2.1 [2|3] 1.2", "C 2 3.1 [3|1] 1.2"]
        @test Gauss._cell_lines(D, 2; bar = true) == ["C 2 [2|3] 1.2", "C 2 [3|1] 1.2"]
        @test sprint(show, D) == "CellTree: 6 cells, by degree [1, 3, 2], LCMOrdering(ComplementDeglex())"
        @test ncells(D) == 6 && ncells(D, 1) == 3 && ncells(D, 3) == 0 && ncells(D, -1) == 0
        @test Gauss.cell_range(D, 2) == 5:6 && Gauss.cell_children(D, 1) == 2:4
        # The cell [3|1]: created by the atom 1 from [3], whose list is (1, 2) (both lcms Δ,
        # a₁ kept as the first in the inherited order).
        x = only(y for y in Gauss.cell_range(D, 2) if Gauss.atom_sequence(D, y) == [3, 1])
        @test Gauss.bar_form(D, x) == [3, 1] && Gauss.cell_weight(D, x) == 4 && Gauss.cell_letter(D, x) == 1
        p = Int(D.parent[x])
        @test Gauss.bar_form(D, p) == [3] && collect(Gauss.admissible_atoms(D, p)) == [1, 2]
        @test collect(Gauss.admissible_atoms(D, 1)) == [1, 2, 3]
        @test check_tree(T) == 0 && check_tree(D) == 0
        # The order of the root acts on the atom sequences of the trivial preorder.
        R = CellTree(dual_dihedral_germ(3), LCMOrdering(TrivialPreorder(); root_order = [3, 2, 1]))
        @test [Gauss.cell_atom(R, y) for y in Gauss.cell_children(R, 1)] == [3, 2, 1]
        @test check_tree(R) == 0 && Gauss._degree_counts(R) == [1, 3, 2]
    end

    @testset "Explicit germs" begin
        for G in (classical_dihedral_germ(5), dual_dihedral_germ(6), signed_interval_germ(4, :D, (2, -4, -3, 1)))
            for (v, (S, _)) in TREE_VARIANTS
                T = CellTree(G, S)
                @test check_tree(T) == 0
                # Euler characteristic of the resolution of ℤ: 0.
                @test sum((-1)^n * ncells(T, n) for n in 0:Gauss._top_degree(T)) == 0
            end
        end
    end

    A(n) = Chevie.coxgroup(:A, n)
    B(n) = Chevie.coxgroup(:B, n)
    cases = [
        "classic_A3" => Chevie.BraidMonoid(A(3)),
        "classic_B3" => Chevie.BraidMonoid(B(3)),
        "dual_A3" => Chevie.DualBraidMonoid(A(3); c = collect(1:3)),
        "dual_B3" => Chevie.DualBraidMonoid(B(3); c = collect(1:3)),
        "dual_D4" => Chevie.DualBraidMonoid(Chevie.coxgroup(:D, 4); c = collect(1:4)),
        "G4" => Chevie.DualBraidMonoid(Chevie.crg(4)),
        "G5" => Chevie.DualBraidMonoid(Chevie.crg(5)),
    ]
    @testset "C lines of J0: $case" for (case, M) in cases
        G = germ(M)
        classical = startswith(case, "classic_")
        for (v, (S, form)) in TREE_VARIANTS
            T = CellTree(G, S)
            @test check_tree(T) == 0
            # Every preorder gives the minimal resolution of a classical braid monoid.
            classical && @test check_classical(T)
            classical && v == "a" && continue
            mine, ref = tree_c_fingerprints(T), reference_c_fingerprints(case, v)
            @test sort(collect(keys(mine))) == sort(collect(keys(ref)))
            k = form === :atoms ? 2 : form === :bar ? 3 : 1
            @test all(haskey(mine, n) && mine[n][1] == ref[n][1] && mine[n][k] == ref[n][k] for n in keys(ref))
        end
    end

    @testset "Lazy germ" begin
        # The same cells for B⁺(A₃) from the lazy germ, whose simples are numbered as met.
        H = Gauss.CoxeterGerm(Chevie.cartan(A(3)))
        L = CellTree(H)
        T = CellTree(germ(Chevie.BraidMonoid(A(3))))
        @test Gauss._degree_counts(L) == [1, 3, 3, 1] && check_classical(L)
        @test all(sort(Gauss._cell_lines(L, n)) == sort(Gauss._cell_lines(T, n)) for n in 0:3)
    end

    @testset "Exports" begin
        for s in (:CellTree, :ncells, :cell_range, :cell_weight, :cell_letter, :cell_atom,
                  :cell_degree, :cell_children, :admissible_atoms, :bar_form, :atom_sequence)
            @test Base.ispublic(Gauss, s) && !Base.isexported(Gauss, s)
        end
    end
end
