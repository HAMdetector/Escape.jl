@testset "EpitopeMap(::String, ::Vector{String}, ::Vector{Int}, ::Vector{Int})" begin
    epitopes = ["SLYNTVATLY", "LYNTVATL"]
    start = [77, 78]
    stop = [86, 85]
    @test EpitopeMap("test", epitopes, start, stop) isa EpitopeMap

    @test_throws ErrorException EpitopeMap("", ["SLYNTVATL", ""], [77], [84]) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("", ["SLYNTVATLY"], [85], [77]) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("", ["SLYNTVATLY"], [-3], [5]) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("", ["SLYNTVATLY"], [77], [85]) isa EpitopeMap
end