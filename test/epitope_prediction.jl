@testset "parse_netmhc(::String)" begin
    testoutput = joinpath(@__DIR__, "data", "netmhc_testoutput.txt")
    df = Escape.parse_netmhc(testoutput, rank_threshold = 100)

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
    @test prediction_a11[1, :affinity] == df[1, :affinity]
end

@testset "epitope_prediction(::AbstractHLAData)" begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_file = joinpath(@__DIR__, "data", "test.fasta")
    data = HLAData("test", fasta_file, hla_types, missing)

    df = Escape.epitope_prediction(data, rank_threshold = 100)

    @test nrow(df) > 0
end

@testset "epitope_feature_matrix(::AbstractHLAData)" begin
    hla_types = [HLAType(parse_allele("A01", "A02", "B05", "B06", "C03", "C07")),
                 HLAType(parse_allele("A01", "A03", "B03", "B05", "C03", "C07")),
                 HLAType(parse_allele("A01", "A04", "B03", "B05", "C03", "C07")),
                 HLAType(parse_allele("A02", "A04", "B03", "B04", "C03", "C07")),
                 HLAType(parse_allele("A02", "A04", "B03", "B04", "C03", "C07"))]

    fasta_file = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_file, hla_types, missing)

    # predicted epitopes at:
    # A02: 5-14, A03: 4-14, C*07: 6
    alleles = Escape.unique_alleles(hla_types)
    expected = zeros(Float64, Escape.sequence_length(hla_data), length(alleles))

    # A01, A02, A03, A04, B03, B04, B05, B06, C03, C07
     expected[4, :] = [0, 0, 1, 0, 0, 0, 0, 0, 0, 0] # A03 has predicted epitope at pos. 6
     expected[5, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 0]
     expected[6, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
     expected[7, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
     expected[8, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
     expected[9, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
    expected[10, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
    expected[11, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
    expected[12, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
    expected[13, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]
    expected[14, :] = [0, 1, 1, 0, 0, 0, 0, 0, 0, 1]

    @test Escape.epitope_feature_matrix(hla_data) == expected
end