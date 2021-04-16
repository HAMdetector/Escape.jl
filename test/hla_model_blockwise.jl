@testset "Escape.run_blockwise(::HLAModel, ::AbstractHLAData)" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]
    r = Escape.replacements(data)

    partitions = Iterators.partition(r, 5)
    res_blockwise = @suppress Escape.run_blockwise(m, data, partitions = partitions,
        iter = 5, warmup = 5, chains = 1)
    res = @suppress Escape.run(m, data, keep_all_parameters = true, iter = 5, warmup = 5, 
        chains = 1)

    for key in keys(res.sf.result[1])
        if startswith(key, "beta_hla")
            @test key in keys(res_blockwise.sf.result[1])
        end
    end
end