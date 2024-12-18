export run_model

function run_model(
    data::AbstractHLAData;
    save_file::String = "",
    stan_kwargs...
)

    if save_file != ""
        save(save_file, Dict("result" => "test"))
        rm(save_file)
    end

    result = Escape.run(
        Escape.HLAModel{4}(), data; stan_kwargs...
    )

    if save_file != ""
        save(save_file, Dict("result" => result))
    end

    return result
end

function run_model(
    model::HLAModel{T}, data::AbstractHLAData;
    save_file::String = "",
    stan_kwargs...
) where T

    if save_file != ""
        save(save_file, Dict("result" => "test"))
        rm(save_file)
    end

    result = Escape.run(
        model, data; stan_kwargs...
    )

    if save_file != ""
        save(save_file, Dict("result" => result))
    end

    return result
end

function Escape.run(
    model::HLAModel{T}, data::AbstractHLAData;
    stan_kwargs... 
) where T

    !ismissing(stan_input(data)) || error("HLAData does not contain a stan_input.")
    
    input = stan_input(data)
    replacements = Escape.replacements(data)
    alleles = sort(unique_alleles(Escape.hla_types(data), 
        allele_depth = allele_depth(data)))

    model_path = joinpath(@__DIR__, "..", "models", "model_" * string(T) * ".stan")
    
    sf = StanInterface.stan(
        model_path, input;
        stan_args = "adapt delta=0.85 algorithm=hmc engine=nuts max_depth=10",
        refresh = 1, iter = 1000, warmup = 300, stan_kwargs...
    )

    return HLAModelResult(model, data, sf, replacements, alleles)
end
