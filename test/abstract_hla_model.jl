@testset "stan_input(::HLAModel, ::AbstractHLAData, ::Replacement)" begin
    struct TestModel <: HLAModel end

    model = TestModel()
    d = HLADataset("Rousseau").data[1]
    replacement = replacements(d)[1]

    @test stan_input(model, d, replacement) isa Dict
end