@testset "HLAAllele(::AbstractString)" begin
    @test HLAAllele("HLA-B*11:10") == HLAAllele("B", "11", "10", missing, missing, missing)
    @test_throws ErrorException HLAAllele("B*11:10")
    @test_throws ErrorException HLAAllele("HLA-A*1110N")
end

@testset "HLAType(::NTuple{6, HLAAllele})" begin
    @test HLAType(parse_allele("A13", "A21", "B51", "B57", "C03", "C07")) isa HLAType
    @test_throws ErrorException HLAType(parse_allele("A13", "B07"))
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

@testset "Base.hash(::HLAAllele)" begin
    @test hash(parse_allele("A*11"), UInt(1)) == hash(parse_allele("A*11"), UInt(1))
    @test hash(parse_allele("A*11")) == hash(parse_allele("A*11"))
    @test hash(parse_allele("B*11")) != hash(parse_allele("B*12"))
    
    allele = HLAAllele("A", "11", missing, missing, missing, missing)
    @test hash(parse_allele("A*11")) == hash(allele)
    
    alleles = [parse_allele("A*11"), parse_allele("A*11"), parse_allele("A*12")]
    unique_alleles = [parse_allele("A*11"), parse_allele("A*12")]
    @test unique(alleles) == unique_alleles
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

@testset "rand(::HLAAllele, ::Symbol)" begin
    @test typeof(rand(HLAAllele, :A)) == HLAAllele
    @test typeof(rand(HLAAllele, :B)) == HLAAllele
    @test typeof(rand(HLAAllele, :C)) == HLAAllele
    @test_throws ErrorException rand(HLAAllele, :Z)

    @test rand(HLAAllele, :A).gene == "A"
    @test rand(HLAAllele, :B).gene == "B"
    @test rand(HLAAllele, :C).gene == "C"
end