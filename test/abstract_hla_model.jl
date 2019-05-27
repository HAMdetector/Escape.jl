@testset "model_path(::HLAModel)" begin
    struct TestModel <: HLAModel end

    @test_throws ErrorException Escape.model_path(TestModel())
    @test Escape.model_path(BernoulliModel()) isa String
end

@testset "stan_input(::HLAModel, ::AbstractHLAData, ::Replacement)" begin
    struct TestModel <: HLAModel end

    model = TestModel()
    d = HLADataset("Rousseau").data[1]
    replacement = replacements(d)[1]

    @test stan_input(model, d, replacement) isa Dict
end