# Tests of the settings of the display (src/complexes/notation.jl, J1.h): construction and
# checks of Notation, copies, the settings of the session and of an IOContext, the names of
# the atoms and their separator. No monoid is built (rule of the layer of complexes).

using GarsideHomology
using Test

const _NT = GarsideHomology.Complexes

@testset "Notation" begin

    @testset "Construction and checks" begin
        N = Notation()
        @test N == _NT._DEFAULT_NOTATION && sprint(show, N) == "Notation()"
        @test N.atom_names === nothing && N.atom_separator === :auto && N.identity == "e"
        @test N.garside_symbol === nothing && N.bar_separator == "|" && N.empty_cell == "e"
        @test N.empty_cell_in_module == "[]" && N.coefficient_style === :auto && N.differential_symbol == "d"
        L = Notation(atom_names = ["a", "b"], coefficient_style = :times)
        @test L.atom_names == ["a", "b"] && L.coefficient_style === :times
        @test sprint(show, L) == "Notation(atom_names = [\"a\", \"b\"], coefficient_style = :times)"
        @test Notation(L; coefficient_style = :auto) == Notation(atom_names = ["a", "b"])
        @test Notation(atom_names = [:x, :y]).atom_names == ["x", "y"]      # converted to strings
        @test Notation(atom_names = :subscripts).atom_names === :subscripts
        f = a -> "x$a"
        @test Notation(atom_names = f).atom_names === f
        @test Notation(atom_separator = ".").atom_separator == "."
        @test Notation(garside_symbol = "Δ").garside_symbol == "Δ"
        @test Notation(coefficient_style = :parentheses).coefficient_style === :parentheses
        @test hash(Notation(atom_names = ["a"])) == hash(Notation(atom_names = ["a"]))
        @test Notation(atom_names = ["a"]) != Notation(atom_names = ["b"])
        for bad in ((atom_names = :letters,), (atom_names = ["a", ""],), (atom_names = ["a", "a"],),
                    (atom_names = 3,), (atom_separator = :none,), (atom_separator = 1,),
                    (identity = "",), (identity = 'e',), (garside_symbol = "",), (bar_separator = "",),
                    (empty_cell = "",), (empty_cell_in_module = "",), (coefficient_style = :dot,),
                    (coefficient_style = "times",), (differential_symbol = "",), (colour = :red,))
            @test_throws ArgumentError Notation(; bad...)
        end
    end

    @testset "Session and IOContext" begin
        try
            @test notation() == Notation()
            N = set_notation!(atom_names = ["a", "b"])
            @test N == notation() == Notation(atom_names = ["a", "b"])
            # The settings not given are kept.
            @test set_notation!(identity = "1") == Notation(atom_names = ["a", "b"], identity = "1")
            @test notation(IOBuffer()) == notation()
            M = Notation(atom_separator = ".")
            @test notation(IOContext(IOBuffer(), :notation => M)) == M
            @test_throws ArgumentError notation(IOContext(IOBuffer(), :notation => 1))
            @test set_notation!() == Notation() == notation()
        finally
            set_notation!()
        end
    end

    @testset "Names of the atoms and separators" begin
        atom_names(N, n) = _NT._atom_names(N, n)
        @test atom_names(Notation(), 3) == ["1", "2", "3"]
        @test atom_names(Notation(atom_names = :subscripts), 12)[[1, 12]] == ["a₁", "a₁₂"]
        @test atom_names(Notation(atom_names = ["a", "b"]), 3) == ["a", "b", "3"]
        @test atom_names(Notation(atom_names = a -> "s" * _NT._subscript_string(a)), 2) == ["s₁", "s₂"]
        @test_throws ArgumentError atom_names(Notation(atom_names = a -> ""), 1)
        @test_throws ArgumentError atom_names(Notation(atom_names = a -> "x"), 2)    # the same name twice
        @test_throws ArgumentError atom_names(Notation(atom_names = ["2"]), 2)       # "2" for the atoms 1 and 2
        separator(N, n) = _NT._atom_separator(N, atom_names(N, n))
        @test separator(Notation(), 9) == "" && separator(Notation(), 10) == "·"
        @test separator(Notation(atom_names = :subscripts), 120) == ""
        @test separator(Notation(atom_names = ["a", "b′", "c"]), 3) == ""
        @test separator(Notation(atom_names = ["ab", "c"]), 2) == "·"
        @test separator(Notation(atom_separator = "."), 3) == "." && separator(Notation(atom_separator = ""), 12) == ""
        @test _NT._is_self_delimiting("s₁₂") && _NT._is_self_delimiting("b′") && _NT._is_self_delimiting("7")
        @test !_NT._is_self_delimiting("₁") && !_NT._is_self_delimiting("") && !_NT._is_self_delimiting("10")
    end

    @testset "Exports" begin
        for s in (:Notation, :set_notation!, :notation)
            @test Base.isexported(GarsideHomology.Complexes, s) && Base.isexported(GarsideHomology, s)
        end
    end
end
