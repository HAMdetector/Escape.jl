export BernoulliModel, BernoulliResult

struct BernoulliModel <: HLAModel
    chains::Int
    iter::Int
end

struct BernoulliResult <: HLAModelResult
    sf::Stanfit
    alleles::Vector{HLAAllele}
end

BernoulliModel(; chains = 4, iter = 2000) = BernoulliModel(chains, iter)

function run(model::BernoulliModel, replacement::Replacement, data::HLAData)
    stan_path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli")
    sf = stan(stan_path, stan_input(model, replacement, data), chains = model.chains,
              iter = model.iter)
    alleles = unique_alleles(data.hla_types)

    return BernoulliResult(sf, alleles)
end

function stan_input(model::BernoulliModel, replacement::Replacement, data::HLAData)
    y = targets(replacement, data)
    m = hla_matrix(data.hla_types)

    stan_input = Dict("y" => collect(skipmissing(y)), "hla_matrix" => m[.!ismissing.(y), :],
                      "n_entries" => length(collect(skipmissing(y))), 
                      "n_alleles" => size(m)[2])
end

function targets(replacement::Replacement, data::HLAData)
    reader = BioSequences.FASTA.Reader(open(data.fasta_file, "r"))
    t = Vector{Union{Missing, Int}}()
    for record in reader
        symbol = Char(BioSequences.FASTA.sequence(record)[replacement.position])
        if symbol ∈ ('X', '-')
            push!(t, missing)
        elseif symbol == replacement.replacement
            push!(t, 1)
        else
            push!(t, 0)
        end
    end

    close(reader)
    return t
end