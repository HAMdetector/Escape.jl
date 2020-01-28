@testset "indices(::Model2Result)" begin
    ds = Escape.HLADataset("Test")
    res = @suppress Escape.run(
        Escape.Model2(), ds.data[1], 
        mincount = 1, iter = 10, warmup = 10, chains = 1
    )

    @test Escape.indices(res) isa Vector{Tuple{T, T}} where T <: Int
    @test Escape.indices(res)[1] == (1, 1)
    @test Escape.indices(res)[end] == (14, 15)
    @test length(Escape.indices(res)) == 14 * 15
end

@testset "thetas(::Model2Result, ::Int, ::Int)" begin
    ds = Escape.HLADataset("Test")
    data = Escape.stan_input(Escape.Model2(), ds.data[1], mincount = 2)

    sf = @suppress Escape.stan(
        joinpath(@__DIR__, "stan", "model_2_gq"), data,
        iter = 10, warmup = 10, chains = 2
    )
    r = Escape.replacements(ds.data[1], mincount = 2)
    alleles = sort(Escape.unique_alleles(ds.data[1].hla_types, depth = 1))
    res = Escape.Model2Result(sf, alleles, r)

    N = sf.data["N"]
    R = length(r)

    for n in 1:N, r in 1:R
        for c in 1:2
            @test length(Escape.thetas(res, r, n)[c]) == 10
            @test all(0 .<= Escape.thetas(res, r, n)[c] .<= 1)
        end
    end
end

@testset "pointwise_loglikelihoods(::Model2Result, ::Int, ::Int)" begin
    ds = Escape.HLADataset("Test")
    data = Escape.stan_input(Escape.Model2(), ds.data[1], mincount = 2)

    sf = @suppress Escape.stan(
        joinpath(@__DIR__, "stan", "model_2_gq"), data, 
        iter = 10, warmup = 10, chains = 2
    )
    r = Escape.replacements(ds.data[1], mincount = 2)
    alleles = sort(Escape.unique_alleles(ds.data[1].hla_types, depth = 1))
    res = Escape.Model2Result(sf, alleles, r)

    N = sf.data["N"]
    R = length(r)

    for n in 1:N, r in 1:R
        for c in 1:2
            expected = sf.result[c]["log_lik.$r.$n"]
            observed = Escape.pointwise_loglikelihoods(res, r, n)[c]
            @test isapprox(observed, expected, rtol = 0.01)
        end
    end

end