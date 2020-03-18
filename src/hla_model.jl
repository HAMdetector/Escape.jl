function Escape.run(
    model::HLAModel{T}, data::AbstractHLAData;
    p0::Real = 5,
    mincount::Int = 10,
    depth::Int = 1,
    stan_kwargs... 
) where T

    input = stan_input(model, data, mincount = mincount, depth = depth)
    input["p0"] = p0

    sf = StanInterface.stan(
        joinpath(@__DIR__, "..", "data", "stan", "model_$T"), input;
        stan_args = "adapt delta=0.97", 
        stan_kwargs...
    )

    return HLAModelResult(model, data, sf)
end

function stan_input(
    model::AbstractHLAModel, data::AbstractHLAData; 
    depth::Int = 1, mincount::Int = 10
)   
    r = replacements(data, mincount = mincount)
    length(r) > 0 || error("no replacements found.")
    ys = ys_maprect(model, data, r)
    xs, X, phy, Z = xs_maprect(model, data, r, depth)

    d = Dict(
        "N" => size(X)[1], 
        "D" => size(X)[2], 
        "R" => length(r), 
        "ys" => ys, 
        "xs" => xs, 
        "X" => X, 
        "phy" => phy,
        "Z" => Z
    )
    
    return d
end

function xs_maprect(
    model::AbstractHLAModel, data::AbstractHLAData, r::Vector{Replacement}, depth::Int
)
    X = Float64.(hla_matrix(data.hla_types; depth = depth))
    for i in 1:size(X)[2]
        X[:, i] .= (X[:, i] .- mean(X[:, i])) ./ std(X[:, i])
    end

    phy = phylogeny_information(model, data, r)
    Z = epitope_information(model, data, r, depth)

    N = size(X)[1]
    D = size(X)[2]
    R = length(r)

    xs = Matrix{Float64}(-1, R, 3 + N*D + D + N)
    fill!(xs, -1)

    for i in 1:R
        t = targets(r[i], data)
        x = combine_input(X[.!ismissing.(t), :]', Z[i, :], phy[i, .!ismissing.(t)])
        xs[i, 1:length(x)] .= x
    end

    return xs, X, phy, Z
end

function ys_maprect(
    model::AbstractHLAModel, data::AbstractHLAData, r::Vector{Replacement}
)   
    R = length(r)
    N = length(data.hla_types)

    ys =  Matrix{Int}(undef, R, N + 1)
    fill!(ys, -1)
    
    for i in 1:R
        t = targets(r[i], data)
        y = combine_input(filter(x -> !ismissing(x), t))
        ys[i, 1:length(y)] .= y
    end

    return ys
end

function epitope_information(
    model::AbstractHLAModel, data::AbstractHLAData, r::Vector{Replacement}, depth::Int
)
    N = length(unique_alleles(data.hla_types, depth = depth))

    if epitope_information(model) == EpitopesIncluded()
        Z = Escape.epitope_feature_matrix(data, depth = depth)[map(x -> x.position, r), :]
    else
        Z = zeros(Int, length(r), N)
    end

    return Z
end

function phylogeny_information(
    model::AbstractHLAModel, data::AbstractHLAData, r::Vector{Replacement}
)
    R = length(r)
    N = length(data.hla_types)
    tree = ismissing(data.tree) ? phylogenetic_tree(data) : data.tree

    if phylogeny_information(model) == PhylogenyIncluded()
        phy = @distributed hcat for i in 1:R
            phylogeny_probabilities(r[i], data, tree)
        end
        phy = phy'
    else
        phy = zeros(Int, R, N)
    end

    return phy
end

function combine_input(x...)
    return [collect(length.(x)); collect(Iterators.flatten(x))]
end