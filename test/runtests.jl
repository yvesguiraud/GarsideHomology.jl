using Gauss
using Test
import Chevie  # loads the package extension GaussChevieExt

@testset "Gauss.jl" begin
    @testset "Package skeleton" begin
        @test isdefined(Gauss, :Complexes)
        @test parentmodule(Gauss.Complexes) === Gauss
    end
    include("test_complexes.jl")
    include("test_germs.jl")
    include("test_table_germ.jl")
    include("test_coxeter_germ.jl")
    include("test_verify.jl")
    include("test_chevie_germs.jl")
    include("test_orderings.jl")
    include("test_snapshots.jl")
    include("test_tree.jl")
    include("test_projection.jl")
    include("test_resolution.jl")
    include("test_homology.jl")
    include("test_pairs.jl")
    include("test_solver_matrices.jl")
    @testset "Chevie extension" begin
        @test Base.get_extension(Gauss, :GaussChevieExt) isa Module
    end
end
