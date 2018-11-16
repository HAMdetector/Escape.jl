@testset "BernoulliModel(::Int, ::Int)" begin
    @test BernoulliModel(4, 2000) isa BernoulliModel
    @test BernoulliModel() isa BernoulliModel
    @test BernoulliModel().iter == 2000
end

@testset "stan_input(::BernoulliModel, ::HLData, ::Replacement)"  begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types)
    replacement = Replacement("test", 2, 'S')

    input = Escape.stan_input(BernoulliModel(), hla_data, replacement)
    
    @test input == Dict("n_entries" => 5, "n_alleles" => 3, 
                        "hla_matrix" => [2 2 2; 2 2 2; 2 2 2; 2 2 2; 2 2 2],
                        "y" => [0, 0, 1, 1, 1])
end

@testset "run(::BernoulliModel, ::HLAData, ::Replacement)" begin
    hla_types = rand(HLAType, 5)
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types)
    replacement = Replacement("test", 2, 'S')
    
    @test Escape.run(BernoulliModel(), hla_data, replacement) isa BernoulliResult
end