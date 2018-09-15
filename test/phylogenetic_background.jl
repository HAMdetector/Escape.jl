@testset "RateMatrix(::Matrix{T}, ::Vector{String})" begin
    @test RateMatrix([0 0; 0 0], ["1", "2"]) isa RateMatrix{Int}
    @test RateMatrix([0.0 1.0; 2.0 3.0], ["0", "1"]) isa RateMatrix{Float64}
    @test RateMatrix(zeros(3, 3), ["a", "b", "c"]) isa RateMatrix
    @test_throws DimensionMismatch RateMatrix([0 0 0; 0 0 0], ["1", "2", "3"])
    @test_throws DimensionMismatch RateMatrix([0 0; 0 0], ["1", "2", "3"])

    @test size(RateMatrix(zeros(2,2), ["1", "2"])) == (2, 2)
    @test size(RateMatrix(zeros(5,5), ["1", "2", "3", "4", "5"])) == (5, 5)
    @test RateMatrix([1 2; 3 4], ["1", "2"])[1, 2] == 2
    @test RateMatrix([1 2; 3 4], ["1", "2"])[2, 1] == 3

    r = RateMatrix([1 2; 3 4], ["1", "2"])
    r[1, 1] = 5
    @test r[1, 1] == 5

    r = RateMatrix([1 2; 3 4], ["A", "B"])
    r["A", "A"] = 2
    @test r[1, 1] == 2
    @test r["A", "A"] == 2
    @test_throws ErrorException r["C", "A"]
end

@testset "transp(::String, ::String, ::RateMatrix, ::Float64)" begin
    m = [-0.886 0.190 0.633 0.063; 
         0.253 -0.696 0.127 0.316; 
         1.266 0.190 -1.519 0.063; 
         0.253 0.949 0.127 -1.329]
    r = RateMatrix(m, ["A", "C", "G", "T"])

    expected = RateMatrix([0.7079 0.0813 0.1835 0.0271;
                           0.1085 0.7377 0.0542 0.0995;
                           0.3670 0.0813 0.5244 0.0271;
                           0.1085 0.2985 0.0542 0.5387], 
                           ["A", "C", "G", "T"])
    
    for i in ["A", "C", "G", "T"], j in ["A", "C", "G", "T"]
        @test abs(Escape.transp(i, j, r, 0.5) - expected[i,j]) < 0.001
    end
end

@testset "L(::PhylogeneticTree, ::RateMatrix)" begin
    tree = PhylogeneticTree("(A:0.3, ((B:0.1, C:0.1):0.1, D:0.2):0.1);")
    annotation = Dict("A" => "A", "B" => "A", "C" => "A", "D" => "A")

    for leaf in leaves(tree)
        name = get_property(tree, leaf, :name)
        set_property!(tree, leaf, :state, annotation[name])
    end

    function hky85(k, A, C, G, T)
        m = [-sum([C, k*G, T]) C k*G T;
             A -sum([A, G, k*T]) G k*T;
             k*A C -sum([k*A, C, T]) T;
             A k*C G -sum([A, k*C, G])]
        μ = -1/sum([A*m[1,1], C*m[2,2], G*m[3,3], T*m[4,4]])

        return RateMatrix(m * μ, ["A", "C", "G", "T"])
    end

    r = hky85(5, 0.4, 0.3, 0.2, 0.1)

    # compares the result to a known value from Huelsenbeck2012
    @test abs(exp(Escape.L(tree, r)) - 0.199465) < 0.0001
end