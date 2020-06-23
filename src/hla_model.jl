function Escape.run(
    model::HLAModel{T}, data::AbstractHLAData;
    p0::Real = 3,
    mincount::Int = 10,
    depth::Int = 1,
    stan_kwargs... 
) where T

    input = stan_input(model, data, depth = depth, mincount = mincount)
    input["p0"] = p0

    sf = StanInterface.stan(
        joinpath(@__DIR__, "..", "data", "stan", "model_$(T)"), input;
        stan_args = "adapt delta=0.90 algorithm=hmc engine=nuts max_depth=12", 
        stan_kwargs...
    )

    return HLAModelResult(model, data, sf)
end

function stan_input(
    model::AbstractHLAModel, data::AbstractHLAData;
    depth::Int = 1, mincount::Int = 10
)
    X = Float64.(hla_matrix(data.hla_types; depth = depth))
    for i in 1:size(X)[2]
        X[:, i] .= (X[:, i] .- mean(X[:, i])) ./ std(X[:, i])
    end

    r = replacements(data, mincount = mincount)
    phy = phylogeny_information(model, data, r)
    Z = epitope_information(model, data, r, depth)
    y = Int[]
    rs = Int[]
    idx = Int[]

    for i in 1:length(r)
        y_i = targets(r[i], data)
        c = .!ismissing.(y_i)

        append!(y, y_i[c])
        append!(rs, repeat([i], sum(c)))
        append!(idx, findall(x -> x == true, c))
    end

    d = Dict(
        "N" => length(y),
        "S" => size(X)[1],
        "D" => size(X)[2],
        "R" => length(r),
        "X" => X,
        "y" => y,
        "rs" => rs,
        "phy" => phy,
        "Z" => Z,
        "idx" => idx
    )

    return d
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
        p = Progress(R, 1)

        phy = Matrix{Float64}(undef, R, N)
        @threads for i in 1:R
            phy[i, :] = phylogeny_probabilities(r[i], data, tree)
            next!(p)
        end
    else
        phy = zeros(Int, R, N)
    end

    return phy
end