@testset "HLAAnalysis(::String, ::HLAModel, ::Vector{HLADataset})" begin
    @test HLAAnalysis("Bernoulli", BernoulliModel(), HLADataset("Rousseau")) isa HLAAnalysis
    @test HLAAnalysis("BernoulliPhylogeny", BernoulliPhylogenyModel(),
        HLADataset("Rousseau")) isa AbstractHLAAnalysis
end

@testset "run(::AbstractHLAAnalysis, ::String)" begin
    fasta_path = joinpath(dirname(@__DIR__), "test",  "data", "test_large.fasta")
    hla_types = rand(HLAType, 15)
    hla_types[5] = HLAType((parse_allele("A11"), parse_allele("A13"),
                            parse_allele("B53"), parse_allele("B56"),
                            missing, missing))
    hla_types[6] = HLAType((missing, missing, missing, missing, missing, missing))
    hla_data_1 = HLAData("protein_1", fasta_path, hla_types, missing)
    hla_data_2 = HLAData("protein_2", fasta_path, rand(HLAType, 15), missing)

    dataset = HLADataset("Test", [hla_data_1, hla_data_2])
    analysis = HLAAnalysis("testrun", BernoulliPhylogenyModel(), dataset)

    testrun_dir = joinpath(dirname(@__DIR__), "test", "data")

    rm(joinpath(testrun_dir, "testrun"), force = true, recursive = true)
    @suppress Escape.run(analysis, testrun_dir)
    @test isdir(testrun_dir)
    result = FileIO.load(joinpath(testrun_dir, "testrun", "analysis_result.jld2"), 
                         "analysis_result")
    @test result isa HLAAnalysisResult
end