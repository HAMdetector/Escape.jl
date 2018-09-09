export PhylogeneticTree

struct PhylogeneticTree
    graph::MetaDiGraph
end

function PhylogeneticTree(newick_string::String)
    graph = MetaDiGraph()
    add_vertex!(graph)
    set_prop!(graph, 1, :name, "root")

    components = newick_components(newick_string)
    for component in components
        add_to_graph!(graph, 1, component)
    end

    return PhylogeneticTree(graph)
end

function add_to_graph!(graph::AbstractMetaGraph, vertex::Int, s::String)
    name = extract_name(s)
    length = extract_length(s)

    add_vertex!(graph)
    set_prop!(graph, nv(graph), :name, name)
    
    add_edge!(graph, vertex, nv(graph))
    set_prop!(graph, Edge(vertex, nv(graph)), :branch_length, length)

    current_vertex = nv(graph)
    if occursin(')', s)
        components = newick_components(s)
        for c in components
            add_to_graph!(graph, current_vertex, c)
        end
    end
end

function extract_name(s::String)
    substring = occursin(')', s) ? s[findlast(isequal(')'), s) + 1 : end] : s
    substring != "" || return missing

    name = String(split(substring, ":")[1])
    return name != "" ? name : missing
end

function extract_length(s::String)
    substring = occursin(')', s) ? s[findlast(isequal(')'), s) + 1 : end] : s
    occursin(':', substring) || return missing

    length = parse(Float64, split(substring, ":")[end])
end

function newick_components(s::String)
    components = Vector{String}()
    component = ""
    
    start = occursin('(', s) ? findfirst(isequal('('), s) + 1 : 1
    stop = occursin(')', s) ? findlast(isequal(')'), s) - 1 : length(s)

    for c in s[start:stop]
        if c == ',' && count(x -> x == '(', component) == count(x -> x == ')', component)
            push!(components, strip(component))
            component = ""
        else
            component *= c
        end
    end

    push!(components, strip(component))
    return components
end

function newick_string(tree::PhylogeneticTree)
    n = neighbors(tree.graph, 1)
    length(n) == 0 && return "();"
    length(n) == 1 && return "(" * newick_string(tree.graph, n[1]) * ");"

    newick = "(" * newick_string(tree.graph, neighbors(tree.graph, 1)[1])
    for n in neighbors(tree.graph, 1)[2:end]
        newick *= ","
        newick *= newick_string(tree.graph, n)
    end
    newick *= ");"
end

function newick_string(graph::AbstractMetaGraph, v::Int)
    name = get_prop(graph, v, :name)
    name = ismissing(name) ? "" : name
    branch_length = get_prop(graph, inneighbors(graph, v)[1], v, :branch_length)
    branch_length = ismissing(branch_length) ? "" : string(branch_length)
    
    if length(neighbors(graph, v)) == 0
        return string(name, branch_length != "" ? ":" : "", branch_length)
    else
        newick = "(" * newick_string(graph, neighbors(graph, v)[1])
        for n in neighbors(graph, v)[2:end]
            newick *= ","
            newick *= newick_string(graph, n)
        end
        newick *= ")" * name
        newick = !ismissing(length) ? newick * ":$branch_length" : newick
    end

    return newick
end