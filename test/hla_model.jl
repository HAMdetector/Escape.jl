@testset "Escape.run(::HLAModel, data::AbstractHLAData)" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]

    @test @suppress Escape.run(m, data, iter = 5, warmup = 5, chains = 1) isa 
        Escape.HLAModelResult
end

# @testset "Escape.reduce_size!(::StanInterface.Stanfit)" begin
#     ds = Escape.HLADataset("Test")
#     m = Escape.HLAModel{4}()
#     data = ds.data[1]
#
#     res = @suppress Escape.run(m, data, iter = 5, warmup = 5, chains = 1)
#     previous_size = Base.summarysize(res.sf)
#
#     @test "aux1_tau.1" in keys(res.sf.result[1])
#     Escape.reduce_size!(res.sf)
#
#     new_size = Base.summarysize(res.sf)
#     @test !("aux1_tau.1" in keys(res.sf.result[1]))
#     @test new_size < previous_size
# end

@testset "Escape.run_model(::HLAModel, ::AbstractHLAData)" begin
    alignment_file_2_digits = joinpath(@__DIR__, "data", 
        "test_large_annotated_2_digits.fasta")
    tree_file = joinpath(@__DIR__, "data", "phylogeny.tree")

    data = @suppress HLAData(
        alignment_file = alignment_file_2_digits, 
        tree_file = tree_file
    )

    @test @suppress Escape.run_model(
        Escape.HLAModel{3}(), 
        data, iter = 5, warmup = 5, chains = 1
    ) isa Escape.HLAModelResult

    save_file = tempname() * ".jld2"
    @test !isfile(save_file)
    @suppress Escape.run_model(
        Escape.HLAModel{3}(), 
        data, iter = 5, warmup = 5, chains = 1, save_file = save_file
    )
    @test isfile(save_file)
    rm(save_file)
end

@testset "Escape.run_model(::AbstractHLAData)" begin
    alignment_file_2_digits = joinpath(@__DIR__, "data", 
    "test_large_annotated_2_digits.fasta")
    tree_file = joinpath(@__DIR__, "data", "phylogeny.tree")

    data = @suppress HLAData(
        alignment_file = alignment_file_2_digits, 
        tree_file = tree_file
    )

    @test @suppress Escape.run_model(
        data, iter = 5, warmup = 5, chains = 1
    ) isa Escape.HLAModelResult

    save_file = tempname() * ".jld2"
    @test !isfile(save_file)
    @suppress Escape.run_model(
        data, iter = 5, warmup = 5, chains = 1, save_file = save_file
    )
    @test isfile(save_file)
    rm(save_file)
end
