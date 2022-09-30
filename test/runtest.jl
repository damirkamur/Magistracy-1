using Test
import FromFile: @from

@testset verbose = true "Magistracy" begin
    include("test_crs_matrix.jl")
end
