@testset "phylogenetic_tree(::String)" begin
    tree = phylogenetic_tree("(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5);")
    names = [get_prop(tree.graph, x, :name) for x in 1:6]
    branch_lengths = [get_prop(tree.graph, e, :branch_length) for
                        e in collect(edges(bfs_tree(tree.graph, 1)))[1:end]]
    @test names == ["root", "A", "B", "E", "C", "D"]
    @test branch_lengths == [0.1, 0.2, 0.5, 0.3, 0.4]
end

@testset "phylogenetic_tree(::AbstractHLAData)" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    hla_types = rand(HLAType, 5)
    hla_data = HLAData("test", fasta_path, hla_types)

    @test phylogenetic_tree(hla_data) isa PhylogeneticTree
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
        @test Escape.newick_string(phylogenetic_tree(n)) == n
    end
end

@testset "leaves(::PhylogeneticTree)" begin
    fasta_path = joinpath(@__DIR__, "data", "phylogeny.fasta")
    hla_data = HLAData("test", fasta_path, rand(HLAType, 6))
    tree = phylogenetic_tree(hla_data)

    @test leaves(tree) == [2,4,6,8,10,11]
end

@testset "isleaf(v::Int, tree::Phylogenetic)" begin
    fasta_path = joinpath(@__DIR__, "data", "phylogeny.fasta")
    hla_data = HLAData("test", fasta_path, rand(HLAType, 6))
    tree =phylogenetic_tree(hla_data)

    @test isleaf(2, tree)
    @test !isleaf(1, tree)
    @test !isleaf(90, tree)
end

@testset "get_property(::PhylogeneticTree, ::Int, ::Symbol)" begin
    fasta_path = joinpath(@__DIR__, "data", "phylogeny.fasta")
    hla_data = HLAData("test", fasta_path, rand(HLAType, 6))
    tree = phylogenetic_tree(hla_data)

    @test get_property(tree, 1, :name) == "root"
    @test get_property(tree, 2, :name) == get_prop(tree.graph, 2, :name)
    @test get_property(tree, 10, :name) == get_prop(tree.graph, 10, :name)
    @test ismissing(get_property(tree, 100, :name))
    @test ismissing(get_property(tree, 1, :nonexisting))
end

@testset "set_property!(::PhylogeneticTree, ::Int, ::Symbol)" begin
    fasta_path = joinpath(@__DIR__, "data", "phylogeny.fasta")
    hla_data = HLAData("test", fasta_path, rand(HLAType, 6))
    tree = phylogenetic_tree(hla_data)

    @test get_property(tree, 1, :name) == "root"
    set_property!(tree, 1, :name, "new name")
    @test get_property(tree, 1, :name) == "new name"
end

@testset "annotate!(::PhylogeneticTree, ::HLAData, ::Replacement)" begin
    fasta_path = joinpath(@__DIR__, "data", "test.fasta")
    r = Replacement("test", 2, 'S')
    hla_data = HLAData("test", fasta_path, rand(HLAType, 5))
    tree = phylogenetic_tree(hla_data)

    annotation = Dict("1" => "0", "2" => "0", "3" => "1", "4" => "1", "5" => "1")
    annotate!(tree, hla_data, r)

    for k in keys(annotation)
        v = filter(x -> get_property(tree, x, :name) == k, leaves(tree))[1]
        @test get_property(tree, v, :state) == annotation[k]
    end
end

@testset "matching(::PhylogeneticTree, ::HLAData)" begin
    fasta_path = joinpath(@__DIR__, "data", "phylogeny.fasta")
    hla_data = HLAData("test", fasta_path, rand(HLAType, 6))
    tree = phylogenetic_tree(hla_data)

    @test Escape.matching(tree, hla_data)

    set_property!(tree, leaves(tree)[2], :name, "test")
    @test Escape.matching(tree, hla_data) isa ErrorException

    short_hla_data = HLAData("short", joinpath(@__DIR__, "data", "test.fasta"),
                             rand(HLAType, 5))
    @test Escape.matching(tree, short_hla_data) isa DimensionMismatch
end