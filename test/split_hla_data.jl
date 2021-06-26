@testset "SplitHLAData" begin
    data = Escape.HLADataset("Test").data[1]
    @test Escape.SplitHLAData(data.name, data.records, data.hla_types, data.tree,
        data.stan_input, 1, 1, data, 1:14, 1) isa Escape.SplitHLAData

    records = Escape.records(joinpath(@__DIR__, "data", "test.fasta"))
    hla_types = rand(HLAType, 5)
    
    split = Escape.SplitHLAData("test", records, hla_types, data.tree, data.stan_input, 
        1, 1, data, 1:14, 1)
    @test typeof(split) <: AbstractHLAData

    wrong_tree = phylogenetic_tree("(X:1.2,5:1.5,(1:0.5,(2:1.6,3:0.7):0.5):0.6);")
    wrong_tree_2 = phylogenetic_tree("(3:0.3,4:0.4)5:0.5);")

    @test_throws ErrorException Escape.SplitHLAData("test", records, hla_types, wrong_tree, 
        data.stan_input, 1, 1, data, 1:14, 1)    
    @test_throws DimensionMismatch Escape.SplitHLAData("test", records, hla_types, 
        wrong_tree_2, data.stan_input, 1, 1, data, 1:14, 1)
end

@testset "split_hla_data(::HLAData, ::Int)" begin
    data = Escape.HLADataset("Test").data[1]
    r = Escape.replacements(data)

    @test Escape.split_hla_data(data, 2) isa Array{Escape.SplitHLAData}
    @test length(Escape.split_hla_data(data, 2)) == 2
    @test Escape.replacements(Escape.split_hla_data(data, 2)[1]) == r[1:7]
end