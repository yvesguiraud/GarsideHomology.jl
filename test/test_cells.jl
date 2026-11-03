# Tests of the cells for the user (src/cells/cell.jl, J1.h): the cells of the bases of a
# resolution, their letters, weight and atoms, their writing with the settings of the
# display, reading them from their bar form or their letters, and the graded basis of a
# resolution.

using GarsideHomology
using GarsideHomology: Cell, CellBasis, GermElement, GradedBasis, basis, cell_tree
using Test
import Chevie

@testset "Cells" begin
    A2 = germ(Chevie.coxgroup(:A, 2))                    # B⁺(A₂): cells e; 1, 2; 2|12
    P = lcm_resolution(A2)
    x = basis(P[2])[1]

    @testset "Letters, weight, atoms" begin
        @test x isa Cell && x.index == 4 && cell_degree(x) == 2 == length(x)
        @test x[1] == A2(2) && x[end] == A2(1, 2) && collect(x) == [A2(2), A2(1, 2)] && first(x) == A2(2)
        @test firstindex(x) == 1 && x[begin] == A2(2)
        @test eltype(x) == GermElement{typeof(A2)} && eltype(typeof(x)) == GermElement{typeof(A2)}
        @test_throws BoundsError x[3]
        @test_throws BoundsError x[0]
        @test cell_weight(x) == A2(1, 2, 1) && atom_sequence(x) == [A2(2), A2(1)]
        e = basis(P[0])[1]
        @test length(e) == 0 && collect(e) == GermElement{typeof(A2)}[] && isone(cell_weight(e))
        @test cell_degree.(x) == 2                       # a scalar for broadcasting
        @test cell_degree.(basis(P[1])) == [1, 1]
    end

    @testset "Equality and order" begin
        @test x == basis(P[2])[1] && hash(x) == hash(basis(P[2])[1])
        y, z = basis(P[1])
        @test y < z < x && sort([x, z, y]) == [y, z, x] && length(Set([y, y, z])) == 2
        Q = lcm_resolution(A2)
        @test basis(Q[2])[1] != x
        @test_throws ArgumentError basis(Q[2])[1] < x
    end

    @testset "Bases of cells" begin
        B = basis(P[1])
        @test B isa CellBasis && B isa AbstractVector{<:Cell} && length(B) == 2 && [c.index for c in B] == [2, 3]
        @test sprint(show, B) == "CellBasis of degree 1 with 2 cells"
        @test sprint(show, MIME"text/plain"(), B) == "2-element CellBasis of degree 1:\n 1\n 2"
        @test B == basis(extend_scalars(P)[1]) && B != basis(lcm_resolution(A2)[1])
        T, U = cell_tree(P), cell_tree(lcm_resolution(A2))
        @test CellBasis(T, 5) == CellBasis(U, 5) && CellBasis(T, 1) != CellBasis(T, 2)   # empty ones are equal
    end

    @testset "Writing" begin
        show_with(N, c) = sprint(show, c; context = :notation => N)
        @test string(x) == "2|12" && string(basis(P[0])[1]) == "e" && repr(basis(P[1])[2]) == "2"
        @test show_with(Notation(atom_names = ["a", "b"]), x) == "b|ab"
        @test show_with(Notation(bar_separator = ", "), x) == "2, 12"
        @test show_with(Notation(empty_cell = "∅"), basis(P[0])[1]) == "∅"
        D65 = lcm_resolution(GarsideHomology._dual_dihedral_germ(65))
        @test string(basis(D65[1])[12]) == "12" && length(basis(D65[1])) == 65
    end

    @testset "Reading" begin
        @test P("2|12") == x == P(" 2 | 12 ") == P("2|1·2") == P(A2(2), A2(1, 2))
        @test P("e") == P() == basis(P[0])[1] && P("1") == basis(P[1])[1]
        @test extend_scalars(P)("2|12") == x                 # the bases are shared
        for bad in ("1|2", "2|3", "", "2||12", "[2|12]", "11")
            @test_throws ArgumentError P(bad)
        end
        @test_throws ArgumentError P(A2(1), A2(2))           # not a cell
        @test_throws ArgumentError P(GarsideHomology.CoxeterGerm([2 -1; -1 2])(1))   # another germ
        @test occursin("between quotes", sprint(showerror, try P(2 | 12) catch e; e end))
        for R in (P, lcm_resolution(GarsideHomology._dual_dihedral_germ(65)),
                  lcm_resolution(germ(Chevie.coxgroup(:B, 3))))
            @test all(R(string(c)) == c for n in range(R) for c in basis(R[n]))
        end
        try
            set_notation!(atom_names = ["a", "b"], bar_separator = " | ")
            @test P("b | ab") == x && string(x) == "b | ab"
        finally
            set_notation!()
        end
    end

    @testset "Graded basis of a resolution" begin
        X = basis(P)
        @test X isa GradedBasis && range(X) == 0:2 && X[2] == basis(P[2]) && X[1][2] == P("2")
        @test sprint(show, X) == "GradedBasis in degrees 0:2, ranks [1, 2, 1]"
        @test sprint(show, MIME"text/plain"(), X) == "GradedBasis in degrees 0:2\n  0: e\n  1: 1, 2\n  2: 2|12"
        @test sprint(show, MIME"text/plain"(), X; context = :notation => Notation(atom_names = ["a", "b"])) ==
              "GradedBasis in degrees 0:2\n  0: e\n  1: a, b\n  2: b|ab"
        @test basis(extend_scalars(P)) == X
        @test Chevie.basis(P) == X && Chevie.basis(P[1]) == basis(P[1])
    end
end
