@testset "Replacement(::String, ::Int, ::Char)" begin
    @test Replacement("Gag", 357, 'S') isa Replacement
end

@testset "targets(::Replacement, ::AbstractHLAData)" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)
    hla_data = HLAData("test", fasta_path, hla_types, missing)

    replacement = Replacement("test", 2, 'S')
    @test all(Escape.targets(replacement, hla_data) .== [0, 0, 1, 1, 1])

    replacement = Replacement("test", 2, 'G')
    @test all(Escape.targets(replacement, hla_data) .== [1, 1, 0, 0, 0])

    replacement = Replacement("test", 3, 'A')
    @test all(Escape.targets(replacement, hla_data) .=== [missing, 1, 1, 1, 1])

    replacement = Replacement("test", 3, '-')
    @test all(Escape.targets(replacement, hla_data) .== [1, 0, 0, 0, 0])

    replacement = Replacement("test", 10, '-')
    @test all(Escape.targets(replacement, hla_data) .== [1, 1, 0, 1, 1])

    replacement = Replacement("test", 10, 'S')
    @test all(Escape.targets(replacement, hla_data) .=== 
        [missing, missing, 1, missing, missing])
end