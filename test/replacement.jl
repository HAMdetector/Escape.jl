@testset "Replacement(::String, ::Int, ::Char, ::Bool)" begin
    @test Replacement("Gag", 357, 'S', true) isa Replacement
end

@testset "targets(::Replacement, ::AbstractHLAData)" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)
    hla_data = HLAData("test", fasta_path, hla_types, missing)

    replacement = Replacement("test", 2, 'S', false)
    @test all(Escape.targets(replacement, hla_data) .== [0, 0, 1, 1, 1])

    replacement = Replacement("test", 2, 'G', false)
    @test all(Escape.targets(replacement, hla_data) .== [1, 1, 0, 0, 0])

    replacement = Replacement("test", 3, 'A', false)
    @test all(Escape.targets(replacement, hla_data) .=== [missing, 1, 1, 1, 1])

    replacement = Replacement("test", 3, '-', false)
    @test all(Escape.targets(replacement, hla_data) .== [1, 0, 0, 0, 0])

    replacement = Replacement("test", 10, '-', false)
    @test all(Escape.targets(replacement, hla_data) .== [1, 1, 0, 1, 1])

    replacement = Replacement("test", 10, 'S', false)
    @test all(Escape.targets(replacement, hla_data) .=== 
        [missing, missing, 1, missing, missing])

    replacement = Replacement("test", 13, 'S', false)
    @test all(Escape.targets(replacement, hla_data) .=== [missing, missing, 1, 1, 1])

    replacement = Replacement("test", 2, 'S', true)
    @test all(Escape.targets(replacement, hla_data) .== [1, 1, 0, 0, 0])
end