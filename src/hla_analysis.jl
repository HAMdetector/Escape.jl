export AbstractHLAAnalysis, HLAAnalysis

abstract type AbstractHLAAnalysis end

struct HLAAnalysis <: AbstractHLAAnalysis
    name::String
    model::HLAModel
    dataset::AbstractHLADataset
end

function run(analysis::AbstractHLAAnalysis, dir::String)
    isdir(dir) || error("$dir must be a directory.")
    to_filename(x::String) = replace(lowercase(x), ' ' => '_')
    
    root = joinpath(dir, to_filename(analysis.name))
    mkdir(root)

    for data in analysis.dataset.data
        data_dir = joinpath(root, to_filename(data.name))
        mkdir(data_dir)

        r = replacements(data)
        @sync @distributed for replacement in r
            result = run(analysis.model, data, replacement)
            filename = string(replacement.protein, "_", replacement.position, 
                replacement.replacement, ".jld2")
            
            FileIO.save(joinpath(data_dir, filename), Dict("result" => result))
        end
    end

    return "done"
end