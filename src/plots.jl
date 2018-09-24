export save, plot

struct GGPlot
    obj::RObject{VecSxp}
end

function save(plot::GGPlot, path::String; width = 14, height = 8.65)
    ggplot_obj = plot.obj
    @rput ggplot_obj path width height

    R"""
    library("ggplot2")
    suppressWarnings(ggsave(path, ggplot_obj, width = width, height = height, unit = "cm"))
    """

    return nothing
end

# WIP
function plot(tree::PhylogeneticTree; title::String = "")
    temp_file = tempname()
    write(temp_file, newick_string(tree))
    n_leaves = length(leaves(tree))
    
    ordered_leaves = [filter(x -> get_property(tree, x, :name) == i, leaves(tree))[1]
                          for i in string.(1:n_leaves)]
    states = [get_property(tree, x, :state) for x in ordered_leaves]
    
    df = state_probabilities(tree)
    probabilities = [df[:prob][findfirst(x -> x == i, df[:sequence])] 
                        for i in string.(1:n_leaves)]

    @rput temp_file n_leaves states title probabilities

    p = R"""
    suppressPackageStartupMessages(library("ggtree"))
    library("tidytree", warn.conflicts = FALSE)

    df = data_frame(label = as.character(1:n_leaves),
                    annotation = states,
                    probabilities = round(probabilities, digits = 2))

    tree <- as.treedata(full_join(as_data_frame(read.newick(temp_file)), 
                                  df, by = "label"))

    p <- ggtree(tree, branch.length = "none", aes(color = annotation)) + 
            #geom_tiplab(aes(label = probabilities), geom = "text") +
            labs(title = title)

    pdf(file = NULL)
    p
    """

    return GGPlot(p)
end