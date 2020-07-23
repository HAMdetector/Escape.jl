@testset "Loo" begin
    if !isfile(joinpath(@__DIR__, "data", "result_loo.jls"))
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test").data[1], 
            mincount = 1, iter = 200, chains = 4
        )

        serialize(joinpath(@__DIR__, "data", "result_loo.jls"), result)
    end

    result = deserialize(joinpath(@__DIR__, "data", "result_loo.jls"))

    sf = Escape.stanfit(result)
    p = extract(sf)

    @test Escape.theta_i(sf, p, 1) isa Vector{Float64}
    @test Escape.pointwise_loglikelihoods(sf, p, 1) isa Vector{Vector{Float64}}
    @test @suppress Escape.loo(result) isa Loo.LooResult
end