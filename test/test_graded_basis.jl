# Tests of the graded basis of a complex, of the generators read from their written form
# and of the printing of long lists (src/complexes/graded_basis.jl, J1.h), on complexes over
# ℤ with bases of descriptors. No monoid is built (rule of the layer of complexes).

using GarsideHomology
using GarsideHomology: basis, GradedBasis
using Test

const _GB = GarsideHomology.Complexes
const _gb_sparse = GarsideHomology.Complexes.SparseArrays.sparse

@testset "Graded bases" begin
    zz = GarsideHomology.Complexes.AbstractAlgebra.ZZ
    M0 = BasedFreeModule(zz, ["p"])
    M1 = BasedFreeModule(zz, ["a", "b"])
    C = ChainComplex([M0, M1], [_gb_sparse([1 -1])])

    @testset "Accessors" begin
        X = basis(C)
        @test X isa GradedBasis && range(X) == 0:1 && length(X) == 2
        @test X[0] == ["p"] && X[1] == ["a", "b"] && X[end] == X[1] && X[firstindex(X)] == X[0]
        @test collect(X) == [["p"], ["a", "b"]] && [n for (n, b) in pairs(X)] == [0, 1]
        @test_throws BoundsError X[2]
        @test X == basis(C) && hash(X) == hash(basis(C))
        @test X != basis(ChainComplex([M0, BasedFreeModule(zz, ["a", "c"])], [_gb_sparse([1 -1])]))
        @test X != basis(ChainComplex([M1], typeof(_gb_sparse([1 -1]))[]))     # other degrees
    end

    @testset "Display" begin
        X = basis(C)
        @test sprint(show, X) == "GradedBasis in degrees 0:1, ranks [1, 2]"
        @test sprint(show, MIME"text/plain"(), X) == "GradedBasis in degrees 0:1\n  0: p\n  1: a, b"
        long = ChainComplex([BasedFreeModule(zz, 65)], typeof(_gb_sparse([1 -1]))[])
        limited = IOContext(IOBuffer(), :limit => true, :displaysize => (24, 40))
        show(limited, MIME"text/plain"(), basis(long))
        @test String(take!(limited.io)) == "GradedBasis in degrees 0:0\n  0: 1, 2, 3, 4, 5, …, 62, 63, 64, 65"
        @test endswith(sprint(show, MIME"text/plain"(), basis(long)), "63, 64, 65")      # not limited: all of them
        line(v, width) = (io = IOContext(IOBuffer(), :limit => true); _GB._print_limited(io, v, width); String(take!(io.io)))
        @test line(1:3, 7) == "1, 2, 3" && line(1:65, 3) == "1, …" && line(["abcdef", "b", "c"], 4) == "abcdef, …"
        @test line(1:5, 13) == "1, 2, 3, 4, 5" && line(1:6, 13) == "1, 2, …, 6"
        @test line(Int[], 10) == ""
    end

    @testset "Generators read from their written form" begin
        @test C("a") == "a" && C(" b ") == "b" && C("p") == "p"
        @test_throws ArgumentError C("z")
        twice = ChainComplex([BasedFreeModule(zz, ["a"]), BasedFreeModule(zz, ["a"])], [_gb_sparse([0;;])])
        @test_throws ArgumentError twice("a")                  # in two degrees
        @test_throws ArgumentError C(14)
        @test occursin("between quotes", sprint(showerror, try C(2 | 12) catch e; e end))
    end
end
