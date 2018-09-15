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
function plot(tree::PhylogeneticTree)
    temp_file = tempname()
    write(temp_file, newick_string(tree))

    @rput temp_file
    p = R"""
    suppressPackageStartupMessages(library("ggtree"))
    library("tidytree", warn.conflicts = FALSE)

    df = data_frame(label = as.character(1:1),
                    annotation = as.character(1:1))
    tree <- as.treedata(full_join(as_data_frame(read.newick(temp_file)), 
                                  df, by = "label"))
    p <- ggtree(tree) +
        geom_tiplab(aes(label = annotation), geom = "text")

    pdf(file = NULL)
    p
    """

    return GGPlot(p)
end