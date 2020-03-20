@testset "RateMatrix(::Matrix{T}, ::Vector{String})" begin
    @test RateMatrix(SA[-1 1; 1 -1], ["1", "2"]) isa RateMatrix
    @test RateMatrix(SA[-1.0 1.0; 1.0 -1.0], ["0", "1"]) isa RateMatrix
    @test RateMatrix(SMatrix{3,3}(zeros(3, 3)), ["a", "b", "c"]) isa RateMatrix
    @test_throws DimensionMismatch RateMatrix(SA[0 0 0; 0 0 0], ["1", "2", "3"])
    @test_throws DimensionMismatch RateMatrix(SA[0 0; 0 0], ["1", "2", "3"])

    @test size(RateMatrix(SA[-1 1; 1 -1], ["1", "2"])) == (2, 2)
    @test size(RateMatrix(SA[-1 0.5 0.5; 0.5 -1 0.5; 0.5 0.5 -1], ["1", "2", "3"])) == (3, 3)
    @test RateMatrix(SA[-2 2; 3 -3], ["1", "2"])[1, 2] == 2
    @test RateMatrix(SA[-4 4; 3 -3], ["1", "2"])[2, 1] == 3

    r = RateMatrix(SA[-1 1; 3 -3], ["A", "B"])
    @test r[1, 1] == -1
    @test r["A", "A"] == -1
    @test_throws ErrorException r["C", "A"]
end

@testset "transp(::RateMatrix, ::String, ::String, ::Float64)" begin
    m = SA[-0.886 0.190 0.633 0.063; 
            0.253 -0.696 0.127 0.316; 
            1.266 0.190 -1.519 0.063; 
            0.253 0.949 0.127 -1.329]
    r = RateMatrix(m, ["A", "C", "G", "T"])

    expected = SA[0.7079 0.0813 0.1835 0.0271;
                  0.1085 0.7377 0.0542 0.0995;
                  0.3670 0.0813 0.5244 0.0271;
                  0.1085 0.2985 0.0542 0.5387]
    
    bases = ["A", "C", "G", "T"]
    for i in bases, j in bases
        base_i_index = findfirst(isequal(i), bases)
        base_j_index = findfirst(isequal(j), bases)
        @test abs(Escape.transp(r, i, j, 0.5) - 
                  expected[base_i_index, base_j_index]) < 0.001
    end
end