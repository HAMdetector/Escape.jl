@testset "indices(::Model4Result)" begin
    ds = Escape.HLADataset("Test")
    res = @suppress Escape.run(
        Escape.Model4(), ds.data[1], 
        mincount = 1, iter = 10, warmup = 10, chains = 2
    )

    @test Escape.indices(res) isa Vector{Tuple{T, T}} where T <: Int
    @test Escape.indices(res)[1] == (1, 1)
    @test Escape.indices(res)[end] == (14, 15)
    @test length(Escape.indices(res)) == 14 * 15
end