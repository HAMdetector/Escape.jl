export PhylogeneticTree, phylogenetic_tree, calculate_tree, newick_string, leaves, isleaf, 
    set_property!, get_property, annotate!

struct PhylogeneticTree
    graph::MetaDiGraph
end

graph(tree::PhylogeneticTree) = getfield(tree, :graph)

function phylogenetic_tree(newick_string::String)
    graph = MetaDiGraph()
    add_vertex!(graph)
    set_prop!(graph, 1, :name, "root")

    components = newick_components(newick_string)
    for component in components
        add_to_graph!(graph, 1, component)
    end

    return PhylogeneticTree(graph)
end

function phylogenetic_tree(data::AbstractHLAData)
    !ismissing(data.phylogenetic_tree) && return data.phylogenetic_tree

    temp_fasta = tempname() * ".fasta"
    temp_name = basename(tempname())
    write_numbered_fasta(data, temp_fasta)

    model = "PROTGAMMAAUTO"
    @suppress Base.run(`raxmlHPC -s $temp_fasta -m $model -p 53 
                        -n $(temp_name) -w $(tempdir()) -T $(Threads.nthreads())`)
    unrooted_tree_path = joinpath(tempdir(), "RAxML_result.$(temp_name)")
    @suppress Base.run(`raxmlHPC -f I -m $model -t $unrooted_tree_path 
                        -w $(tempdir()) -n $(temp_name * "_rooted") 
                        -T $(Threads.nthreads())`)

    final_tree_path = joinpath(tempdir(), "RAxML_rootedTree.$(temp_name * "_rooted")")
    
    return phylogenetic_tree(readline(final_tree_path))
end

function write_numbered_fasta(data::AbstractHLAData, filepath::String)
    records = Escape.records(data)

    write_numbered_fasta(records, filepath)
end

function write_numbered_fasta(records::Vector{FASTX.FASTA.Record}, filepath::String)
    open(FASTA.Writer, filepath) do writer
        for (i, record) in enumerate(records)
            identifier = string(i)
            sequence = FASTA.sequence(record)
            new_record = FASTA.Record(identifier, sequence)

            write(writer, new_record)
        end
    end
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
    graph = Escape.graph(tree)
    n = LightGraphs.neighbors(graph, 1)
    length(n) == 0 && return "();"
    length(n) == 1 && return "(" * newick_string(graph, n[1]) * ");"

    newick = "(" * newick_string(graph, LightGraphs.neighbors(graph, 1)[1])
    for n in LightGraphs.neighbors(graph, 1)[2:end]
        newick *= ","
        newick *= newick_string(graph, n)
    end
    newick *= ");"

    return newick
end

function newick_string(graph::AbstractMetaGraph, v::Int)
    name = get_prop(graph, v, :name)
    name = ismissing(name) ? "" : name
    branch_length = get_prop(graph, inneighbors(graph, v)[1], v, :branch_length)
    branch_length = ismissing(branch_length) ? "" : string(branch_length)
    
    if length(LightGraphs.neighbors(graph, v)) == 0
        return string(name, branch_length != "" ? ":" : "", branch_length)
    else
        newick = "(" * newick_string(graph, LightGraphs.neighbors(graph, v)[1])
        for n in LightGraphs.neighbors(graph, v)[2:end]
            newick *= ","
            newick *= newick_string(graph, n)
        end
        newick *= ")" * name
        newick = !ismissing(length) ? newick * ":$branch_length" : newick
    end

    return newick
end

function isleaf(v::Int, tree::PhylogeneticTree)
    v ∈ leaves(tree) && return true
    return false
end

function leaves(tree::PhylogeneticTree)
    graph = Escape.graph(tree)
    leaves = map(x -> x[1], attracting_components(graph))

    return leaves
end

function get_property(tree::PhylogeneticTree, v::Int, property::Symbol)
    graph = Escape.graph(tree)
    has_prop(graph, v, property) || return missing

    return get_prop(graph, v, property)
end

function set_property!(tree::PhylogeneticTree, v::Int, property::Symbol, value)
    graph = Escape.graph(tree)

    set_prop!(graph, v, property, value)
end

function annotate!(tree::PhylogeneticTree, data::AbstractHLAData, replacement::Replacement)
    name(data) == protein(replacement) || error("Replacement does not match HLA data.")
    records = Escape.records(data)

    for (i, record) in enumerate(records)
        sequence = FASTA.sequence(record)
        symbol = Char(sequence[position(replacement)])

        v = filter(x -> get_property(tree, x, :name) == string(i), leaves(tree))[1]
        
        if !negated(replacement)
            if symbol == Escape.replacement(replacement)
                set_property!(tree, v, :state, "1")
            elseif symbol in ('-', 'X')
                set_property!(tree, v, :state, missing)
            else
                set_property!(tree, v, :state, "0")
            end
        elseif negated(replacement)
            if symbol == Escape.replacement(replacement)
                set_property!(tree, v, :state, "0")
            elseif symbol in ('-', 'X')
                set_property!(tree, v, :state, missing)
            else
                set_property!(tree, v, :state, "1")
            end
        end
    end
end

function matching(tree::PhylogeneticTree, data::AbstractHLAData)
    records = Escape.records(data)
    matching = Escape.matching(tree, records)

    return matching
end

function matching(tree::PhylogeneticTree, records::Vector{FASTX.FASTA.Record})
    n_leaves = length(leaves(tree))

    msg = string(
        "Number of leaves ($n_leaves) does not match number of sequences ",
        "($(length(records)))"
    )
    n_leaves == length(records) || return DimensionMismatch(msg)

    msg = """:name property of leaves must be "1", "2", ..., "$n_leaves"."""
    leaf_names = Set(map(x -> get_property(tree, x, :name), leaves(tree)))
    leaf_names == Set(string.(1:n_leaves)) || return ErrorException(msg)

    return true
end

function states(tree::PhylogeneticTree)
    states = Vector{String}(undef, length(leaves(tree)))

    for leaf in leaves(tree)
        state = get_property(tree, leaf, :state)
        name = get_property(tree, leaf, :name)

        states[parse(Int, name)] = ismissing(state) ? "-" : state
    end

    return states
end