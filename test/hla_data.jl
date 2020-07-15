@testset "HLAData(::String, ::String, ::Vector{HLAType})" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)

    @test typeof(HLAData("test", fasta_path, hla_types, missing, missing)) <: AbstractHLAData
    @test HLAData("test", fasta_path, hla_types, missing, missing) isa HLAData

    tree = phylogenetic_tree("(4:1.2,5:1.5,(1:0.5,(2:1.6,3:0.7):0.5):0.6);")
    @test typeof(HLAData("test", fasta_path, hla_types, tree, missing)) <: AbstractHLAData
    @test HLAData("test", fasta_path, hla_types, tree, missing) isa HLAData

    wrong_tree = phylogenetic_tree("(X:1.2,5:1.5,(1:0.5,(2:1.6,3:0.7):0.5):0.6);")
    wrong_tree_2 = phylogenetic_tree("(3:0.3,4:0.4)5:0.5);")

    @test_throws ErrorException HLAData("test", fasta_path, hla_types, wrong_tree, missing) 
    @test_throws DimensionMismatch HLAData("test", fasta_path, hla_types, wrong_tree_2, missing)
end