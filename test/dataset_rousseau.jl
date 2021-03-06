@testset "HLADataset(Rousseau)" begin
    @test HLADataset("Rousseau") isa HLADataset

    dataset = HLADataset("Rousseau")
    @test dataset.name == "Rousseau"
    @test dataset.data[1] isa HLAData
    @test dataset.data[2] isa HLAData
    @test dataset isa AbstractHLADataset
end