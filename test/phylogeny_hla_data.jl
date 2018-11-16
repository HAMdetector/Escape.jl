@testset "PhylogenyHLAData(::String, ::String, ::Vector{HLAData}, ::PhylogeneticTree)" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)
    tree = phylogenetic_tree("(1:0.1,2:0.2,(3:0.3,4:0.4)5:0.5);")
    hla_data = PhylogenyHLAData("test", fasta_path, hla_types, tree)

    @test hla_data isa PhylogenyHLAData
    @test hla_data isa AbstractHLAData
end

@testset "phylogenetic_tree(::PhylogenyHLAData)" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)
    tree = phylogenetic_tree("(1:0.1,2:0.2,(3:0.3,4:0.4)5:0.5);")
    hla_data = PhylogenyHLAData("test", fasta_path, hla_types, tree)

    @test phylogenetic_tree(hla_data) == tree
end