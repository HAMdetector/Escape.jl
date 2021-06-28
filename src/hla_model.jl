function run_model(
    data::AbstractHLAData;
    keep_all_parameters::Bool = false,
    save_file::String = "",
    stan_kwargs...
) where T

    if save_file != ""
        save(save_file, Dict("result" => "test"))
        rm(save_file)
    end

    result = Escape.run(
        Escape.HLAModel{4}(), data; keep_all_parameters = keep_all_parameters, stan_kwargs...
    )

    if save_file != ""
        save(save_file, Dict("result" => result))
    end

    return result
end

function run_model(
    model::HLAModel{T}, data::AbstractHLAData;
    keep_all_parameters::Bool = false,
    save_file::String = "",
    stan_kwargs...
) where T

    if save_file != ""
        save(save_file, Dict("result" => "test"))
        rm(save_file)
    end

    result = Escape.run(
        model, data; 
        keep_all_parameters = keep_all_parameters, stan_kwargs...
    )

    if save_file != ""
        save(save_file, Dict("result" => result))
    end

    return result
end

function Escape.run(
    model::HLAModel{T}, data::AbstractHLAData;
    keep_all_parameters::Bool = false,
    stan_kwargs... 
) where T

    !ismissing(stan_input(data)) || error("HLAData does not contain a stan_input.")
    
    input = stan_input(data)
    replacements = Escape.replacements(data)
    alleles = sort(unique_alleles(Escape.hla_types(data), 
        allele_depth = allele_depth(data)))

    sf = StanInterface.stan(
        joinpath(@__DIR__, "..", "models", "model_$(T)"), input;
        stan_args = "adapt delta=0.85 algorithm=hmc engine=nuts max_depth=10",
        refresh = 1, iter = 1000, warmup = 300, stan_kwargs...
    )

    !keep_all_parameters && reduce_size!(sf)

    return HLAModelResult(model, data, sf, replacements, alleles)
end

function reduce_size!(sf::StanInterface.Stanfit)
    for d in sf.result
        for key in keys(d)
            if !any(startswith.(key, ["lp__", "b_phy", "b0_hla", "beta_hla"]))
                delete!(d, key)
            end
        end
    end

    return sf
end