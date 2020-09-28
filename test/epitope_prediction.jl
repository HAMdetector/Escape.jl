@testset "sequence_length(::Vector{FASTX.FASTA.Record})" begin
    fasta_file = joinpath(@__DIR__, "data", "test.fasta")
    records = FASTX.FASTA.Record[]

    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            push!(records, record)
        end
    end

    @test Escape.sequence_length(records) == 14
end

@testset "consensus_sequence(::Vector{FASTX.FASTA.Record})" begin
    fasta_file = joinpath(@__DIR__, "data", "test.fasta")
    records = FASTX.FASTA.Record[]

    open(FASTA.Reader, fasta_file) do reader
        for record in reader
            push!(records, record)
        end
    end

    @test Escape.consensus_sequence(records) == "MSARASVMGSRASV"
end

@testset "epitope_prediction(::String)" begin
    peptide = "ARGDEFFPE"
    alleles = Escape.valid_alleles()
    prediction = Escape.epitope_prediction(peptide, alleles[1:10])

    @test prediction isa DataFrame
end

@testset "epitope_prediction(::AbstractHLAData)" begin
    hla_types = [HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01")),
                 HLAType(parse_allele("A01", "A01", "B01", "B01", "C01", "C01"))]
    fasta_file = joinpath(@__DIR__, "data", "test.fasta")
    data = HLAData("test", fasta_file, hla_types, missing, missing)

    df = @suppress Escape.epitope_prediction(data, rank_threshold = 100)

    @test nrow(df) > 0
end

@testset "epitope_feature_matrix(::AbstractHLAData)" begin
    hla_types = [HLAType(parse_allele("A01", "A02", "B05", "B06", "C03", "C07")),
                 HLAType(parse_allele("A01", "A03", "B03", "B05", "C03", "C07")),
                 HLAType(parse_allele("A01", "A04", "B03", "B05", "C03", "C07")),
                 HLAType(parse_allele("A02", "A04", "B03", "B04", "C03", "C07")),
                 HLAType(parse_allele("A02", "A04", "B03", "B04", "C03", "C07"))]

    fasta_file = joinpath(@__DIR__, "data", "test.fasta")
    hla_data = HLAData("test", fasta_file, hla_types, missing, missing)

    @test Escape.epitope_feature_matrix(hla_data) isa Matrix
end