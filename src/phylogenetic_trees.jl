export PhylogeneticTree

struct PhylogeneticTree
    graph::MetaDiGraph
end

function PhylogeneticTree(newick_string::String)
    graph = MetaDiGraph()
    add_vertex!(graph)
    set_prop!(graph, 1, :name, "root")
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