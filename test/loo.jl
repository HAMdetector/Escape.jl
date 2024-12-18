@testset "Loo" begin
    result = @suppress Escape.run(
        Escape.HLAModel{4}(), Escape.HLADataset("Test").data[1], 
        iter = 10, warmup = 10, chains = 4, seed = 123
    )

    sf = Escape.stanfit(result)
    p = StanInterface.extract(sf)

    @test Escape.theta_i(sf, p, 1) isa Vector{Float64}
    @test Escape.pointwise_loglikelihoods(sf, p, 1) isa Vector{Vector{Float64}}
    @test @suppress Escape.loo(result) isa Loo.LooResult
end
