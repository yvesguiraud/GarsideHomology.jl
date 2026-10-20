using Gauss
using Test
import Chevie  # loads the package extension GaussChevieExt

@testset "Gauss.jl" begin
    @testset "Package skeleton" begin
        @test isdefined(Gauss, :Complexes)
        @test parentmodule(Gauss.Complexes) === Gauss
    end
    include("test_complexes.jl")
    @testset "Chevie extension" begin
        @test Base.get_extension(Gauss, :GaussChevieExt) isa Module
    end
end
