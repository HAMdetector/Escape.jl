@testset "Loo" begin
    tmp = tempname()
    try
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test"),
            mincount = 1,
            result_dir = tempdir(),
            result_name = splitdir(tmp)[end],
            iter = 200, warmup = 200, chains = 4
        )

        sf = Escape.stanfit(result[1])
        p = extract(sf)

        @test Escape.theta_i(sf, p, 1) isa Vector{Float64}
        @test Escape.pointwise_loglikelihoods(sf, p, 1) isa Vector{Vector{Float64}}
        @test @suppress Escape.loo(result[1]) isa Loo.LooResult
        @test @suppress Escape.loo(result) isa Loo.LooResult
    finally
        rm(tmp, force = true, recursive = true)
    end
end