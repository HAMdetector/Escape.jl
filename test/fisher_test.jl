@testset "FisherTest()" begin
    @test FisherTest() isa FisherTest
end

@testset "Escape.run(::FisherTest, ::AbstractHLAData, ::Replacement)" begin
    hla_types = [HLAType(parse_allele("A01", "A02", "B05", "B06", "C03", "C07")),
                 HLAType(parse_allele("A01", "A03", "B03", "B05", "C03", "C07")),
                 HLAType(parse_allele("A01", "A04", "B03", "B05", "C03", "C07")),
                 HLAType(parse_allele("A02", "A04", "B03", "B04", "C03", "C07")),
                 HLAType(parse_allele("A02", "A04", "B03", "B04", "C03", "C07"))]

    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    replacement = Replacement("test", 2, 'S', false)
    hla_data = HLAData("test", fasta_path, hla_types, missing)

    y = Bool[0, 0, 1, 1, 1] # representation for replacement
    x = Bool[1, 1, 1, 0, 0] # allele vector for A01

    result = Escape.run(FisherTest(), hla_data, replacement)
    counts, p, beta = Escape.fisher_exact_test(y, x)

    @test result.counts[parse_allele("A01")] == counts
    @test result.p_values[parse_allele("A01")] == p
    @test result.log_odds[parse_allele("A01")] == beta 
end

@testset "fisher_exact_test(::Vector{Bool}, ::Vector{Bool}" begin
    a = Bool[1, 1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 1, 1] # y's
    b = Bool[0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0] # alleles

    #             with allele  without allele           
    # with y    |      1              7
    # without y |      4              2                            
    expected_result = [1 7; 4 2]
    counts, p, beta = Escape.fisher_exact_test(a, b)

    @test counts == [1 7; 4 2]
    @test p == pvalue(FisherExactTest(1, 7, 4, 2), method = :minlike)
    @test beta == log(FisherExactTest(1, 7, 4, 2).ω)
end