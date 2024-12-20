@testset "accessor functions for HLAModelResult" begin
    tmp = tempname()

    if !isfile(joinpath(@__DIR__, "data", "result.jls"))
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test").data[1], 
            iter = 10, chains = 1
        )

        serialize(joinpath(@__DIR__, "data", "result.jls"), result)
    end

    result = deserialize(joinpath(@__DIR__, "data", "result.jls"))

    @test Escape.stan_input(result) == StanInterface.stan_data(Escape.stanfit(result))
    @test Escape.stanfit(result) === result.sf
    @test Escape.hla_data(result) === result.data
    @test Escape.hla_model(result) === result.model
end

@testset "replacement_summary(::HLAModelResult)" begin
    if !isfile(joinpath(@__DIR__, "data", "result.jls"))
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test").data[1],
            result_dir = tempdir(),
            result_name = splitdir(tmp)[end],
            iter = 200, warmup = 200, chains = 4
        )

        serialize(joinpath(@__DIR__, "data", "result.jls"), result)
    end

    result = deserialize(joinpath(@__DIR__, "data", "result.jls"))

    @test @suppress Escape.replacement_summary(result) isa DataFrame
end

# @testset "reduce_size!(::HLAModelResult)" begin
#     ds = Escape.HLADataset("Test")
#     m = Escape.HLAModel{4}()
#     data = ds.data[1]
#
#     res = @suppress Escape.run(m, data, iter = 5, warmup = 5, chains = 1, 
#         keep_all_parameters = true)
#
#     previous_size = Base.summarysize(res)
#     Escape.reduce_size!(res)
#
#     @test Base.summarysize(res) < previous_size
# end

@testset "diagnostics(::HLAModelResult)" begin
    ds = Escape.HLADataset("Test")
    data = ds.data[1]

    res = @suppress Escape.run(Escape.HLAModel{1}(), data, iter = 5, warmup = 5, chains = 1)

    @test @suppress Escape.diagnostics(res) == StanInterface.diagnose(res.sf)
end

@testset "load_result(::String)" begin
    ds = Escape.HLADataset("Test")

    save_file = tempname() * ".jld2"
    res = @suppress Escape.run_model(ds.data[1], warmup = 10, iter = 10, chains = 1, 
        save_file = save_file)

    @test Escape.load_result(save_file) isa Escape.HLAModelResult
end
