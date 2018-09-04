@testset "HLAData(::String, ::String, ::Vector{HLAType})" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)

    @test typeof(HLAData("test", fasta_path, hla_types)) <: AbstractHLAData
    @test HLAData("test", fasta_path, hla_types) isa HLAData
end

@testset "hla_matrix(::Vector{NTuple{6, HLAAllele}})" begin
    hla_types = [HLAType(parse_allele("A11", "A03", "B7", "B7", "C7", "C7")),
                 HLAType(parse_allele("A33", "A03", "B7", "B7", "C7", "C7"))]

    @test Escape.hla_matrix(hla_types) == [1 1 0 2 2; 1 0 1 2 2]
end