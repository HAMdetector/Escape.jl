@testset "Escape.loo(::HLAModelResult)" begin
    ds = Escape.HLADataset("Test")
    res = @suppress Escape.run(
        Escape.Model2(), ds.data[1], 
        mincount = 1, iter = 10, warmup = 10, chains = 2
    )

    @test Escape.loo(res) isa Loo.LooResult
end