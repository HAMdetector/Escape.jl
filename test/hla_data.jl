@testset "HLAData(::String, ::String, ::Vector{HLAType})" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)

    @test typeof(HLAData("test", fasta_path, hla_types, missing, missing, 
        allele_depth = 1, replacement_mincount = 1)) <: AbstractHLAData
    @test HLAData("test", fasta_path, hla_types, missing, missing, 
        allele_depth = 1, replacement_mincount = 1) isa HLAData

    tree = phylogenetic_tree("(4:1.2,5:1.5,(1:0.5,(2:1.6,3:0.7):0.5):0.6);")
    @test typeof(HLAData("test", fasta_path, hla_types, tree, missing, 
        allele_depth = 1, replacement_mincount = 1)) <: AbstractHLAData
    @test HLAData("test", fasta_path, hla_types, tree, missing, 
        allele_depth = 1, replacement_mincount = 1) isa HLAData

    wrong_tree = phylogenetic_tree("(X:1.2,5:1.5,(1:0.5,(2:1.6,3:0.7):0.5):0.6);")
    wrong_tree_2 = phylogenetic_tree("(3:0.3,4:0.4)5:0.5);")

    @test_throws ErrorException HLAData(
        "test", fasta_path, hla_types, wrong_tree, missing, 
        allele_depth = 1, replacement_mincount = 1
    ) 
    @test_throws DimensionMismatch HLAData(
        "test", fasta_path, hla_types, wrong_tree_2, missing, 
        allele_depth = 1, replacement_mincount = 1
    )

    @test_throws ErrorException HLAData(
        "test", fasta_path, hla_types, tree, missing, 
        allele_depth = 3, replacement_mincount = 1
    )
end

@testset "check_hla_matrix(::AbstractMatrix)" begin
    @test_throws ErrorException Escape.check_hla_matrix(zeros(50, 0))
    @test_logs (:warn,) Escape.check_hla_matrix([0 0 1; 0 0 0; 1 1 1])
    @test isnothing(Escape.check_hla_matrix([0 0 1; 1 1 1; 1 0 1]))
end

@testset "check_inputs(::String, ::String))" begin
    alignment_file = joinpath(@__DIR__, "data", "test_large_annotated_2_digits.fasta")
    tree_file = joinpath(@__DIR__, "data", "phylogeny.tree")

    @test isnothing(Escape.check_inputs(alignment_file, tree_file))
    @test_throws ArgumentError Escape.check_inputs(
        joinpath(@__DIR__, "data", "test.fasta"), tree_file
    )

    @test_throws ArgumentError Escape.check_inputs(
        alignment_file, joinpath(@__DIR__, "data", "test.tree")
    )
end

@testset "check_inputs(::String, ::String, ::String)" begin
    alignment_file = joinpath(@__DIR__, "data", "test_large.fasta")
    tree_file = joinpath(@__DIR__, "data", "phylogeny.tree")
    hla_annotation_file = joinpath(@__DIR__, "data", "hla_annotation_2_digits.csv")

    @test_throws ArgumentError Escape.check_inputs(
        joinpath(@__DIR__, "data", "no_alignment.fasta"), tree_file, hla_annotation_file
    )
    @test_throws ArgumentError Escape.check_inputs(
        alignment_file, joinpath(@__DIR__, "data", "no_phylogeny.tree"), hla_annotation_file
    )
    @test_throws ArgumentError Escape.check_inputs(
        alignment_file, tree_file, joinpath(@__DIR__, "data", "no_annotation.csv")
    )
    @test_throws ArgumentError Escape.check_inputs(
        joinpath(@__DIR__, "data", "test.fasta"), tree_file, hla_annotation_file
    )
    @test_throws ArgumentError Escape.check_inputs(
        joinpath(@__DIR__, "data", "test_large_wrong_id.fasta"), 
        tree_file, 
        hla_annotation_file
    )
    @test_throws ArgumentError Escape.check_inputs(
        alignment_file, 
        joinpath(@__DIR__, "data", "phylogeny_id_wrong_leaf.tree"), 
        hla_annotation_file
    )

    @test Escape.check_inputs(alignment_file, tree_file, hla_annotation_file) |> isnothing

    tree_file = joinpath(@__DIR__, "data", "phylogeny_id_leaves.tree")
    @test Escape.check_inputs(alignment_file, tree_file, hla_annotation_file) |> isnothing
end

@testset "hla_annotation_df(::String)" begin
    
    @test_throws ErrorException Escape.hla_annotation_df(
        joinpath(@__DIR__, "data", "hla_annotation_wrong_format.csv")
    )
    @test_logs (:warn,) Escape.hla_annotation_df(
        joinpath(@__DIR__, "data", "hla_annotation_sparse.csv")
    )
    @test Escape.hla_annotation_df(
        joinpath(@__DIR__, "data", "hla_annotation_2_digits.csv")
    ) isa DataFrame
    @test Escape.hla_annotation_df(
        joinpath(@__DIR__, "data", "hla_annotation_4_digits.csv")
    ) isa DataFrame
end

@testset "phylogeny_information(::AbstractHLAData, ::Vector{Replacement})" begin
    ds = Escape.HLADataset("Test")
    data = ds.data[1]
    r = Escape.replacements(ds.data[1])

    @test @suppress collect(Escape.phylogeny_information(data, r)) isa Matrix
end

@testset "epitope_information(::AbstractHLAData, ::Vector{Replacement}, ::Int)" begin
    ds = Escape.HLADataset("Test")
    data = ds.data[1]
    r = Escape.replacements(ds.data[1])

    @test @suppress collect(Escape.epitope_information(data, r, 1)) isa Matrix
end

@testset "compute_stan_input(::AbstractHLAData)" begin
    ds = Escape.HLADataset("Test")
    data = ds.data[1]

    @test @suppress Escape.compute_stan_input(data) isa Dict
end

@testset "HLAData(; kwargs...)" begin
    alignment_file_2_digits = joinpath(@__DIR__, "data", 
        "test_large_annotated_2_digits.fasta")
    alignment_file_4_digits = joinpath(@__DIR__, "data", 
        "test_large_annotated_4_digits.fasta")
    alignment_file_ID = joinpath(@__DIR__, "data", "test_large.fasta")
    tree_file = joinpath(@__DIR__, "data", "phylogeny.tree")
    hla_annotation_file_2_digits = joinpath(@__DIR__, "data", "hla_annotation_2_digits.csv")
    hla_annotation_file_4_digits = joinpath(@__DIR__, "data", "hla_annotation_4_digits.csv")


    @test @suppress HLAData(
        alignment_file = alignment_file_2_digits, 
        tree_file = tree_file,
        replacement_mincount = 2
    ) isa Escape.HLAData

    @test @suppress HLAData(
        alignment_file = alignment_file_4_digits, 
        tree_file = tree_file,
        allele_depth = 2
    ) isa Escape.HLAData

    @test @suppress HLAData(
        alignment_file = alignment_file_ID, 
        tree_file = tree_file, 
        hla_annotation_file = hla_annotation_file_2_digits,
    ) isa Escape.HLAData

    @test @suppress HLAData(
        alignment_file = alignment_file_ID,
        tree_file = tree_file,
        hla_annotation_file = hla_annotation_file_4_digits,
        allele_depth = 1
    ) isa Escape.HLAData

    @test @suppress HLAData(
        alignment_file = alignment_file_ID,
        tree_file = tree_file,
        hla_annotation_file = hla_annotation_file_4_digits,
        allele_depth = 2
    ) isa Escape.HLAData
end