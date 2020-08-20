function Escape.run(
    model::HLAModel{T}, data::AbstractHLAData;
    p0::Real = 10,
    mincount::Int = 1,
    depth::Int = 1,
    stan_kwargs... 
) where T

    input = stan_input(model, data, depth = depth, mincount = mincount)
    replacements = Escape.replacements(data, mincount = mincount)
    alleles = sort(unique_alleles(Escape.hla_types(data), depth = depth))
    input["p0"] = p0

    sf = StanInterface.stan(
        joinpath(@__DIR__, "..", "data", "stan", "model_$(T)"), input;
        stan_args = "adapt delta=0.85 algorithm=hmc engine=nuts max_depth=10",
        stan_kwargs...
    )

    return HLAModelResult(model, data, sf, replacements, alleles)
end

function stan_input(
    model::AbstractHLAModel, data::AbstractHLAData;
    depth::Int = 1, mincount::Int = 1
)
    ismissing(data.stan_input) || return data.stan_input

    hla_types = Escape.hla_types(data)
    X = Float64.(hla_matrix_standardized(hla_types; depth = depth))
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
        "idx" => idx,
        "p0" => 10
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
    tree = phylogenetic_tree(data)

    if phylogeny_information(model) == PhylogenyIncluded()
        phy = @showprogress "phylogeny " @distributed hcat for i in 1:R
            phylogeny_probabilities(r[i], data, tree)
        end
        phy = phy'
    else
        phy = zeros(Int, R, N)
    end

    return phy
end