@testset "HLAAnalysis(::String, ::HLAModel, ::Vector{HLADataset})" begin
    @test HLAAnalysis("Bernoulli", BernoulliModel(), HLADataset("Rousseau")) isa HLAAnalysis
    @test HLAAnalysis("BernoulliPhylogeny", BernoulliPhylogenyModel(),
        HLADataset("Rousseau")) isa AbstractHLAAnalysis
end