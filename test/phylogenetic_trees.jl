@testset "PhylogeneticTree(::String)" begin
    tree = PhylogeneticTree("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);")
    names = [get_prop(tree.graph, x, :name) for x in 1:6]
    branch_lengths = [get_prop(tree.graph, e, :branch_length) for
                        e in collect(edges(bfs_tree(tree.graph, 1)))[1:end]]
    @test names == ["root", "A", "B", "E", "C", "D"]
    @test branch_lengths == [0.1, 0.2, 0.5, 0.3, 0.4]
end

@testset "add_to_graph!(::AbstractMetaGraph, ::Int, ::String)" begin
    graph = MetaDiGraph()
    add_vertex!(graph)

    Escape.add_to_graph!(graph, 1, "(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)")
    names = [get_prop(graph, x, :name) for x in [3,4,5,6,7]]
    branch_lengths = [get_prop(graph, e, :branch_length) for 
                        e in collect(edges(graph))[2:end]]
    
    @test names == ["A", "B", "E", "C", "D"]
    @test branch_lengths == [0.1, 0.2, 0.5, 0.3, 0.4]

end

@testset "extract_name(::String)" begin
    @test ismissing(Escape.extract_name(""))
    @test ismissing(Escape.extract_name("(A:3, B:1)"))
    @test ismissing(Escape.extract_name("(A:3, B:1):3"))
    @test ismissing(Escape.extract_name("(,)"))
    @test Escape.extract_name("B:5") == "B"
    @test Escape.extract_name("(A:3, B:1)E") == "E"
    @test Escape.extract_name("(A:3, B:1)E:5") == "E"
    @test Escape.extract_name("(,)E") == "E"
    @test Escape.extract_name("(,B)E:3") == "E"
end

@testset "extract_length(::String)" begin
    @test ismissing(Escape.extract_length(""))
    @test ismissing(Escape.extract_length("(A:3, B:1)"))
    @test ismissing(Escape.extract_length("(,)"))
    @test ismissing(Escape.extract_length("(A:3, B:1)E"))
    @test ismissing(Escape.extract_length("(,)E"))
    @test Escape.extract_length("(A:3, B:1):3") == 3.0
    @test Escape.extract_length("B:5") == 5.0
    @test Escape.extract_length("(A:3, B:1)E:5") == 5.0
    @test Escape.extract_length("(,B)E:3") == 3.0
end

@testset "newick_component(::String)" begin
    components = ["A:0.1", "B:0.2", "(C:0.3,D:0.4)E:0.5"]
    @test Escape.newick_components("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);") == components
    @test Escape.newick_components("(,)") == ["", ""]
    @test Escape.newick_components("(,B)E:3") == ["", "B"]
end

@testset "newick_string(::PhylogeneticTree)" begin
    newick = ["(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);", "(,);",
              "((B:0.2,(C:0.3,D:0.4)E:0.5)A:0.1);", 
              "(B:6.0,(A:5.0,C:3.0,E:4.0):5.0,D:11.0);"]
    for n in newick
        @test Escape.newick_string(PhylogeneticTree(n)) == n
    end
end