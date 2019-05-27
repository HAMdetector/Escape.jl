struct TestModel <: HLAModel end
struct UndefinedPath <: HLAModel end
struct TestPhylogenyModel <: HLAPhylogenyModel end

@testset "stan_input(::HLAPhylogenyModel, ::AbstractHLAData, ::Replacement)" begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    r = Replacement("test", 2, 'S')

    d = HLAData("test", fasta_path, hla_types, missing)
    tree = phylogenetic_tree(d)

    without_tree = d
    with_tree = HLAData("test", fasta_path, hla_types, tree)

    @test stan_input(TestPhylogenyModel(), without_tree, r) isa Dict
    @test stan_input(TestPhylogenyModel(), with_tree, r) isa Dict
end

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