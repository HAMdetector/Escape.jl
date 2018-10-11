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