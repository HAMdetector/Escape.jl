@testset "accessor functions for HLAModelResult" begin
    ds = HLADataset("Test")
    res = @suppress Escape.run(Escape.HLAModel{4}(), ds.data[1], mincount = 1, 
        iter = 10, chains = 1)

    @test Escape.stan_input(res) === res.sf.data
    @test Escape.stanfit(res) === res.sf
    @test Escape.hla_data(res) === res.data
    @test Escape.hla_model(res) === res.model
end

@testset "replacement_summary(::HLAModelResult)" begin
    ds = HLADataset("Test")
    res = @suppress Escape.run(Escape.HLAModel{4}(), ds.data[1], mincount = 1, 
        iter = 10, chains = 1)

    @test Escape.replacement_summary(res) isa DataFrame
end