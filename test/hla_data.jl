@testset "HLAData(::String, ::String, ::Vector{HLAType})" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)

    @test typeof(HLAData("test", fasta_path, hla_types)) <: AbstractHLAData
    @test HLAData("test", fasta_path, hla_types) isa HLAData
end