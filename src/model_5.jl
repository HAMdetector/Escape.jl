struct Model5 <: HLAModel end

struct Model5Result <: HLAModelResult 
    sf::Stanfit
    alleles::Vector{HLAAllele}
    replacements::Vector{Replacement}
end

function hallo()
    return "bla"
end

function stan_input(
    model::Model5, data::AbstractHLAData; 
    depth::Int = 1, mincount::Int = 10
)   
    println("in")
    tree = ismissing(data.tree) ? phylogenetic_tree(data) : data.tree
    X = hla_matrix(data.hla_types; depth = depth)
    r = replacements(data, mincount = mincount)
    Z = Escape.epitope_feature_matrix(data)[map(x -> x.position, r), :]
    N = size(X)[1]
    D = size(X)[2]
    R = length(r)

    phy = Vector{Float64}(undef, length(r))
    for i in 1:R
        phy[i] = phylogeny_probabilities(r[i], data, tree)
    end

    # for i in 1:R
    #     t = targets(r[i], data)
    #     y_mean[i] = mean(skipmissing(t))
    #     count = 0
    #     for j in 1:N
    #         if !ismissing(t[j])
    #             count += 1
    #             ys[i, count + 1] = t[j]
    #         end
    #     end

    #     ys[i, 1] = count
    #     X_ = hcat(phy[.!ismissing.(t), i], X[.!ismissing.(t), :])
    #     xs[i, 1:((D + 1) * count)] = reshape(X_', :, 1)
    # end


    # d = Dict("N" => N, "D" => D, "R" => R, "ys" => ys, "xs" => xs,
    #          "y_mean" => y_mean, "hla_mean" => hla_mean, "p0" => 5,
    #          "Z" => Z)
    
    # return d
end