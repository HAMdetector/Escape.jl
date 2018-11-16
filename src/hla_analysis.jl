export AbstractHLAAnalysis, HLAAnalysis, AbstractHLAAnalysisResult, HLAAnalysisResult

abstract type AbstractHLAAnalysis end
abstract type AbstractHLAAnalysisResult end

struct HLAAnalysis <: AbstractHLAAnalysis
    name::String
    model::HLAModel
    dataset::AbstractHLADataset
end

struct HLAAnalysisResult <: AbstractHLAAnalysisResult
    analysis::HLAAnalysis
    path::String
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
        if analysis.model isa HLAPhylogenyModel
            tree = PhylogeneticTree(data)
        end

        @sync @distributed for replacement in r
            if analysis.model isa HLAPhylogenyModel
                result = run(analysis.model, data, replacement, tree)
            else
                result = run(analysis.model, data, replacement)
            end

            filename = string(replacement.protein, "_", replacement.position, 
                replacement.replacement, ".jld2")
            
            FileIO.save(joinpath(data_dir, filename), Dict("result" => result))
        end
    end

    analysis_result = HLAAnalysisResult(analysis, root)
    FileIO.save(joinpath(analysis_result.path, "analysis_result.jld2"),
                Dict("analysis_result" => analysis_result))
end