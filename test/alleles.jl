@testset "HLAAllele(::AbstractString)" begin
    @test HLAAllele("HLA-B*11:10") == HLAAllele("B", "11", "10", missing, missing, missing)
    @test_throws ErrorException HLAAllele("B*11:10")
    @test_throws ErrorException HLAAllele("HLA-A*1110N")
end

@testset "HLAType(::NTuple{6, HLAAllele})" begin
    @test HLAType(parse_allele("A13", "A21", "B51", "B57", "C03", "C07")) isa HLAType
    @test_throws MethodError HLAType(parse_allele("A13", "B07"))
end

@testset "hla_matrix(::Vector{HLAType})" begin
    hla_types = [HLAType(parse_allele("A11", "A03", "B7", "B7", "C7", "C7")),
                 HLAType(parse_allele("A33", "A03", "B7", "B7", "C7", "C7"))]

    @test Escape.hla_matrix(hla_types) == [1 1 0 1 1; 1 0 1 1 1]
end

@testset "hla_matrix_standardized(::Vector{HLAType})" begin
    hla_types = [HLAType(parse_allele("A11", "A03", "B7", "B7", "C7", "C7")),
                 HLAType(parse_allele("A33", "A03", "B7", "B7", "C7", "C7"))]

    @test Escape.hla_matrix_standardized(hla_types) ≈ [0 0.5 -0.5 0 0; 0 -0.5 0.5 0 0]
end

@testset "Base.in(::HLAAllele, ::HLAType)" begin
    hla_type = HLAType(parse_allele("A13", "A21:01", "B51", "B57", "C03", "C07"))
    @test parse_allele("A13") in hla_type
    @test parse_allele("A21:01") in hla_type
    @test parse_allele("B51") ∈ hla_type
    @test parse_allele("B52") ∉ hla_type
    @test parse_allele("A21") in hla_type
end

@testset "Base.:(==)(::HLAAllele, ::HLAAllele)" begin
    @test HLAAllele("HLA-B*11:10") == HLAAllele("HLA-B*11:10")
    @test HLAAllele("HLA-B*11:10") != HLAAllele("HLA-B*11:11")
end

@testset "is_valid_allele(::HLAAllele)" begin
    @test is_valid_allele("HLA-B") == false
    @test is_valid_allele("A*11:01") == false
    @test is_valid_allele("HLA-A*11:01") == true
    @test is_valid_allele("HLA-A*11:01N") == true
end

@testset "Base.isless(::HLAAllele, ::HLAAllele)" begin
    @test parse_allele("A11") > parse_allele("A03")
    @test parse_allele("A11") < parse_allele("B11")
    @test parse_allele("A11") < parse_allele("C20")
    @test parse_allele("A03") <= parse_allele("A11")
    @test parse_allele("A11") <= parse_allele("A11")
    @test ismissing(parse_allele("A11:01") <= parse_allele("A11"))
    @test parse_allele("B51:03") > parse_allele("B51:01")
end

@testset "parse_allele(::AbstractString)" begin
    @test parse_allele("A2402") == HLAAllele("A", "24", "02", missing, missing, missing)
    @test parse_allele("A2401") != parse_allele("A2402")
    @test parse_allele("Cw*0302") == HLAAllele("C", "03", "02", missing, missing, missing)
    @test parse_allele("Cw0302") == parse_allele("HLA-C*03:02")
    @test parse_allele("Cw4") == HLAAllele("C", "04", missing, missing, missing, missing)
    @test parse_allele("A*29") == HLAAllele("A", "29", missing, missing, missing, missing)
    @test parse_allele("A021010102N") == HLAAllele("A", "02", "101", "01", "02", "N")
    @test parse_allele("A02(1010102N)") == HLAAllele("A", "02", "101", "01", "02", "N")
    @test parse_allele("A02:10N") == HLAAllele("A", "02", "10", missing, missing, "N")
    @test parse_allele("HLA-B*44:02:01:02S") == HLAAllele("B", "44", "02", "01", "02", "S")
    @test ismissing(parse_allele("?"))
    @test ismissing(parse_allele("N"))
end

@testset "parse_allele(x...)" begin
    @test length(parse_allele("A*11", "A*12")) == 2
    @test parse_allele("A*11", "A*12") == (parse_allele("A*11"), parse_allele("A*12"))
end

@testset "string(::HLAAllele)" begin
    @test string(HLAAllele("HLA-A*11:01:05N")) == "HLA-A*11:01:05N"
    @test string(parse_allele("A1")) == "HLA-A*01"
end

@testset "limit_hla_accuracy(::HLAAllele)" begin
    allele = parse_allele("A*11:01:05N")
    @test limit_hla_accuracy(allele) == parse_allele("A*11")
    @test limit_hla_accuracy(allele, depth = 2) == parse_allele("A*1101")
    @test limit_hla_accuracy(parse_allele("A11")) == parse_allele("A*11")
end

@testset "rand(::HLAAllele, ::Symbol)" begin
    @test typeof(rand(HLAAllele, :A)) == HLAAllele
    @test typeof(rand(HLAAllele, :B)) == HLAAllele
    @test typeof(rand(HLAAllele, :C)) == HLAAllele
    @test_throws ErrorException rand(HLAAllele, :Z)

    @test rand(HLAAllele, :A).gene == "A"
    @test rand(HLAAllele, :B).gene == "B"
    @test rand(HLAAllele, :C).gene == "C"
end

@testset "rand(::HLAType)" begin
    @test rand(HLAType) isa HLAType
    @test length(rand(HLAType, 3)) == 3
    @test rand(HLAType, 3)[2] isa HLAType
end

@testset "unique_alleles" begin
    hla_types = [HLAType(parse_allele("A11", "A03", "B03", "B51", "C03", "C03")),
                 HLAType(parse_allele("A13", "A27", "B03", "B51", "C07", "C07"))]
    @test length(Escape.unique_alleles(hla_types)) == 8
    @test allunique(Escape.unique_alleles(hla_types))
end

@testset "hla_accuracy(::HLAAllele)" begin
    @test hla_accuracy(parse_allele("A11")) == 1
    @test hla_accuracy(parse_allele("A11:01")) == 2
    @test hla_accuracy(parse_allele("B57:01:01")) == 3
    @test hla_accuracy(parse_allele("B57:01:01:01")) == 4
    @test hla_accuracy(parse_allele("B57:01:01:01N")) == 5
end
