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

function run(model::BernoulliModel, data::HLAData, replacement::Replacement)
    path = joinpath(@__DIR__, "..", "data", "stan", "bernoulli")
    input = stan_input(model, data, replacement)
    sf = stan(path, input, chains = model.chains, iter = model.iter)
    alleles = sort(unique_alleles(filter(x -> missing ∉ x, data.hla_types)))

    return BernoulliResult(sf, alleles)
end

function stan_input(model::BernoulliModel, data::HLAData, replacement::Replacement)
    y = targets(data, replacement)
    
    complete_cases = findall(i -> !ismissing(y[i]) && missing ∉ data.hla_types[i], 
        1:length(y))
    y = y[complete_cases]
    m = hla_matrix(data.hla_types[complete_cases])

    stan_input = Dict("y" => y, "hla_matrix" => m,
                      "n_entries" => length(y), 
                      "n_alleles" => size(m)[2])
end

function targets(data::HLAData, replacement::Replacement)
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