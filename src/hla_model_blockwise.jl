function run_partition(model::HLAModel{T}, data::AbstractHLAData;
    mincount::Int = 1,
    depth::Int = 1,
    partitions = Iterators.partition(replacements(data, mincount = mincount), 100),
    wp::AbstractWorkerPool = default_worker_pool(),
    stan_kwargs...
) where T

    input = stan_input(model, data, depth = depth, mincount = mincount)
    alleles = sort(unique_alleles(hla_types(data), depth = depth))
    
    sf_files = pmap(x -> run_partition_(model, data, x; 
        mincount = mincount, depth = depth, wp = wp, stan_kwargs...), wp, partitions)
    sf_1 = deserialize(sf_files[1]) # used to access iter, chains field etc.
    combined_d = [Dict{String, Vector{Float64}}() for i in 1:sf_1.chains]

    i = 0
    for file in sf_files
        sf = deserialize(file)
        for r_idx in 1:sf.data["R"]
            i = i + 1
            for d in 1:sf.data["D"]
                for (res_idx, res) in enumerate(sf.result)
                    combined_d[res_idx]["b_phy.$i"] = res["b_phy.$r_idx"]
                    combined_d[res_idx]["b0_hla.$i"] = res["b0_hla.$r_idx"]
                    combined_d[res_idx]["beta_hla.$i.$d"] = res["beta_hla.$r_idx.$d"]
                end
            end
        end
    end

    combined_sf = StanInterface.Stanfit(sf_1.model, sf_1.data, sf_1.iter, sf_1.chains,
        combined_d, sf_1.diagnostics)

    rm.(sf_files)

    return HLAModelResult(model, data, combined_sf, vcat(partitions...), alleles)
end

function run_partition_(model::HLAModel{T}, data::AbstractHLAData, replacements;
    mincount, depth, wp, stan_kwargs...
) where T

    input = stan_input(model, data, depth = depth)
    alleles = sort(unique_alleles(hla_types(data), depth = depth))
    all_replacements = Escape.replacements(data, mincount = mincount)
    
    idx_(replacement) = findfirst(x -> x == replacement, all_replacements)
    r_idx = map(idx_, replacements)
    
    y = Int[]
    rs = Int[]
    idx = Int[]

    for i in 1:length(replacements)
        y_i = targets(replacements[i], data)
        c = .!ismissing.(y_i)

        append!(y, y_i[c])
        append!(rs, repeat([i], sum(c)))
        append!(idx, findall(x -> x == true, c))
    end

    d = Dict(
        "N" => length(y),
        "S" => input["S"],
        "D" => input["D"],
        "R" => length(replacements),
        "X" => input["X"],
        "y" => y,
        "rs" => rs,
        "phy" => input["phy"][r_idx, :],
        "Z" => input["Z"][r_idx, :],
        "idx" => idx
    )

    filename, io = mktemp()
    sf = StanInterface.stan(
        joinpath(@__DIR__, "..", "models", "model_$(T)"), d;
        stan_args = "adapt delta=0.85 algorithm=hmc engine=nuts max_depth=10",
        refresh = 1, stan_kwargs...
    )

    serialize(io, sf)
    close(io)

    return filename
end