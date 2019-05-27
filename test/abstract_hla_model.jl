struct TestModel <: HLAModel end
struct UndefinedPath <: HLAModel end

@testset "run(::HLAModel, ::AbstractHLAData, ::Replacement)" begin
    Escape.model_path(::TestModel) = 
        joinpath(@__DIR__, "..", "data", "stan", "bernoulli_hs")

    d = HLADataset("Rousseau").data[1]
    r = replacements(d)[1]

    @test Escape.model_path(TestModel()) == 
        joinpath(@__DIR__, "..", "data", "stan", "bernoulli_hs")

    result = @suppress Escape.run(TestModel(), d, r, chains = 1, iter = 100)
    @test result isa Escape.GenericHLAModelResult
end

@testset "model_path(::HLAModel)" begin
    @test_throws ErrorException Escape.model_path(UndefinedPath())
    @test Escape.model_path(BernoulliModel()) isa String
end

@testset "stan_input(::HLAModel, ::AbstractHLAData, ::Replacement)" begin
    model = TestModel()
    d = HLADataset("Rousseau").data[1]
    replacement = replacements(d)[1]

    @test stan_input(model, d, replacement) isa Dict
end