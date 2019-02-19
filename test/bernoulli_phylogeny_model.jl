@testset "BernoulliPhylogenyModel(::Int, ::Int)" begin
    @test BernoulliPhylogenyModel(2000, 4) isa BernoulliPhylogenyModel
    @test BernoulliPhylogenyModel() isa BernoulliPhylogenyModel
    @test BernoulliPhylogenyModel().iter == 2000
end

@testset "stan_input(::BernoulliModel, ::HLAData, ::Replacement, ::PhylogeneticTree)"  begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")

    hla_data = HLAData("test", fasta_path, hla_types, missing)
    tree = phylogenetic_tree(hla_data)
    replacement = Replacement("test", 2, 'S')
    annotate!(tree, hla_data, replacement)

    input = Escape.stan_input(BernoulliPhylogenyModel(), hla_data, replacement, tree)
    p = Escape.state_probabilities(tree, TwoStateGTR)
    phylogeny_effect = [p[s]["1"] for s in string.(1:5)]

    @test input == Dict("n_entries" => 5, "n_alleles" => 3, 
                        "hla_matrix" => [2 2 2; 2 2 2; 2 2 2; 2 2 2; 2 2 2],
                        "y" => [0, 0, 1, 1, 1],
                        "phylogeny_effect" => phylogeny_effect)
end

@testset "run(::BernoulliPhylogenyModel, ::HLAData, ::Replacement)" begin
    hla_types = rand(HLAType, 5)
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types, missing)
    replacement = Replacement("test", 2, 'S')
    
    @test @suppress Escape.run(BernoulliPhylogenyModel(), hla_data, replacement) isa 
        BernoulliPhylogenyResult
end

@testset "BernoulliPhylogenyModel with finnish horseshoe prior" begin
    hla_types = rand(HLAType, 5)
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types, missing)
    hla_data.tree = phylogenetic_tree(hla_data)

    replacement = Replacement("test", 2, 'S')

    model = BernoulliPhylogenyModel(prior = :finnish_horseshoe)
    result = @suppress Escape.run(model, hla_data, replacement) 
    @test result isa BernoulliPhylogenyResult
end

@testset "BernoulliPhylogenyModel with broad t prior" begin
    hla_types = rand(HLAType, 5)
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types, missing)
    hla_data.tree = phylogenetic_tree(hla_data)
    replacement = Replacement("test", 2, 'S')

    model = BernoulliPhylogenyModel(prior = :broad_t)
    result = @suppress Escape.run(model, hla_data, replacement) 
    @test result isa BernoulliPhylogenyResult
end