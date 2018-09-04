@testset "Replacement(::String, ::Int, ::Char)" begin
    @test Replacement("Gag", 357, 'S') isa Replacement
end