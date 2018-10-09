@testset "BernoulliPhylogenyModel(::Int, ::Int)" begin
    @test BernoulliPhylogenyModel(4, 2000) isa BernoulliPhylogenyModel
    @test BernoulliPhylogenyModel() isa BernoulliPhylogenyModel
    @test BernoulliPhylogenyModel().iter == 2000
end

@testset "stan_input(::BernoulliModel, ::PhylogeneticTree, ::Replacement, ::HLData)"  begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types)
    tree = PhylogeneticTree(hla_data)
    replacement = Replacement("test", 2, 'S')
    annotate!(tree, hla_data, replacement)

    input = Escape.stan_input(BernoulliPhylogenyModel(), tree, replacement, hla_data)
    p = Escape.state_probabilities(tree, TwoStateGTR)
    phylogeny_effect = [p[s]["1"] for s in string.(1:5)]

    @test input == Dict("n_entries" => 5, "n_alleles" => 3, 
                        "hla_matrix" => [2 2 2; 2 2 2; 2 2 2; 2 2 2; 2 2 2],
                        "y" => [0, 0, 1, 1, 1],
                        "phylogeny_effect" => phylogeny_effect)
end

@testset "run(::BernoulliModel, ::Replacement, ::HLAData)" begin
    hla_types = rand(HLAType, 5)
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types)
    replacement = Replacement("test", 2, 'S')
    
    @test @suppress Escape.run(BernoulliPhylogenyModel(), replacement, hla_data) isa 
        BernoulliPhylogenyResult
end

@testset "run(::BernoulliModel, ::Phylogenetictree, ::Replacement, ::HLAData)" begin
    hla_types = rand(HLAType, 5)
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types)
    tree = PhylogeneticTree(hla_data)
    replacement = Replacement("test", 2, 'S')
    
    @test @suppress Escape.run(BernoulliPhylogenyModel(), tree, replacement, hla_data) isa 
        BernoulliPhylogenyResult
end