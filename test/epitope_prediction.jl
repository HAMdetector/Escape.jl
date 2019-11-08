@testset "parse_netmhc(::String)" begin
    testoutput = joinpath(@__DIR__, "data", "netmhc_testoutput.txt")
    df = Escape.parse_netmhc(testoutput)

    @test df isa DataFrame
    @test nrow(df) == 1
end

@testset "epitope_prediction(::String)" begin
    testoutput = joinpath(@__DIR__, "data", "netmhc_testoutput.txt")
    df = Escape.parse_netmhc(testoutput)

    peptide = "ARGDEFFPE"
    prediction = Escape.epitope_prediction(peptide)

    @test prediction isa DataFrame
    @test nrow(prediction) > 0

    prediction_a11 = filter(x -> x[:allele] == parse_allele("A1101"), prediction)
    @test prediction_a11[!, :affinity] == df[!, :affinity]
end

@testset "epitope_prediction_fasta(::String)" begin
    fastafile = joinpath(@__DIR__, "data", "test.fasta")
    df = Escape.epitope_prediction_fasta(fastafile)

    @test df isa DataFrame
    @test nrow(df) > 0
end

@testset "epitope_prediction_fasta(::String, ::Int)" begin
    fastafile = joinpath(@__DIR__, "data", "test.fasta")
    df = Escape.epitope_prediction_fasta(fastafile, 1)

    @test df isa DataFrame
    @test nrow(df) > 0
end

@testset "epitope_prediction(::Replacement, ::HLAData)" begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_path, hla_types, missing)
    replacement = Replacement("test", 2, 'S')

    df = Escape.epitope_prediction(replacement, hla_data)
    
    @test df isa DataFrame
    @test nrow(df) > 0
end