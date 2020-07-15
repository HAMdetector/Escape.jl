@testset "EpitopeMap(::String, ::Vector{String}, ::Vector{Int}, ::Vector{Int})" begin
    epitopes = ["SLYNTVATLY", "LYNTVATL"]
    start = [77, 78]
    stop = [86, 85]
    @test EpitopeMap("test", epitopes, start, stop, 
        [parse_allele("B51"), parse_allele("A01")], 100) isa EpitopeMap

    @test_throws ErrorException EpitopeMap("", ["SLYNTVATL", ""], [77], [84], 
        [parse_allele("A11")], 100) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("", ["SLYNTVATLY"], [85], [77],
        [parse_allele("A03")], 100) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("", ["SLYNTVATLY"], [-3], [5],
        [parse_allele("A01")], 100) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("", ["SLYNTVATLY"], [77], [85],
        [parse_allele("A68")], 100) isa EpitopeMap
    @test_throws ErrorException EpitopeMap("test", epitopes, start, stop, 
        [parse_allele("B51")], 40) isa EpitopeMap
end

@testset "epitope_map(::AbstractHLAData)" begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types, missing, missing)

    map = Escape.epitope_map(hla_data)

    @test map isa Escape.EpitopeMap
    @test map.maplength == 14
end

@testset "Base.in(::Pair{Replacement, HLAAllele}, ::EpitopeMap)" begin
    epitopes = ["SLYNTVATLY", "LYNTVATL"]
    start = [77, 78]
    stop = [86, 85]
    alleles = [parse_allele("B51"), parse_allele("A01")]
    map = EpitopeMap("test", epitopes, start, stop, alleles, 100)
    r_1 = Replacement("test", 77, 'G', false)
    r_2 = Replacement("test", 76, 'S', false)
    r_3 = Replacement("test", 99, 'S', false)
    r_4 = Replacement("test", 80, 'G', false)

    @test (r_1 => parse_allele("B51")) in map
    @test !((r_1 => parse_allele("A01")) in map)
    @test !((r_1 => parse_allele("B52")) in map)
    @test !((r_2 => parse_allele("B51")) in map)
    @test !((r_2 => parse_allele("B51")) in map)
    @test !((r_3 => parse_allele("B51")) in map)
    @test !((r_3 => parse_allele("B52")) in map)
    @test ((r_4 => parse_allele("B51")) in map)
    @test ((r_4 => parse_allele("A01")) in map)
    @test !((r_4 => parse_allele("B52")) in map)

end