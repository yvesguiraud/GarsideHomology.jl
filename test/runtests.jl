using GarsideHomology
using Test
import Chevie  # loads the package extension GarsideHomologyChevieExt

# verbose: the summary shows the time of each file of tests, and of the parts of the
# longest one, the germs from Chevie (J1.g, task g3).
@testset "GarsideHomology.jl" verbose = true begin
    @testset "Package skeleton" begin
        @test isdefined(GarsideHomology, :Complexes)
        @test parentmodule(GarsideHomology.Complexes) === GarsideHomology
    end
    include("test_complexes.jl")
    include("test_notation.jl")
    include("test_graded_basis.jl")
    include("test_module_elements.jl")
    include("test_germs.jl")
    include("test_table_germ.jl")
    include("test_coxeter_germ.jl")
    include("test_germ_elements.jl")
    include("test_verify.jl")
    include("test_germ_isomorphism.jl")
    include("test_chevie_germs.jl")
    include("test_orderings.jl")
    include("test_snapshots.jl")
    include("test_tree.jl")
    include("test_cells.jl")
    include("test_projection.jl")
    include("test_resolution.jl")
    include("test_resolution_elements.jl")
    include("test_certificates.jl")
    include("test_modules.jl")
    include("test_convenience.jl")
    include("test_homology.jl")
    include("test_expected_values.jl")
    include("test_pairs.jl")
    include("test_solver_matrices.jl")
    @testset "Chevie extension" begin
        @test Base.get_extension(GarsideHomology, :GarsideHomologyChevieExt) isa Module
    end
end
