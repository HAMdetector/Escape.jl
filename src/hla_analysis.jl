export AbstractHLAAnalysis, HLAAnalysis, AbstractHLAAnalysisResult, HLAAnalysisResult,
       analysis_result

abstract type AbstractHLAAnalysis end
abstract type AbstractHLAAnalysisResult end

struct HLAAnalysis <: AbstractHLAAnalysis
    name::String
    model::HLAModel
    dataset::AbstractHLADataset
end

struct HLAAnalysisResult <: AbstractHLAAnalysisResult
    analysis::AbstractHLAAnalysis
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

        @sync @distributed for replacement in r
            if analysis.model isa HLAPhylogenyModel
                result = run(analysis.model, data, replacement, phylogenetic_tree(data);
                             wp = WorkerPool([myid()]))
            else
                result = run(analysis.model, data, replacement; wp = WorkerPool([myid()]))
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

function analysis_result(dirpath::AbstractString)
    isdir(dirpath) || error("directory $filepath does not exist.")
    result_path = joinpath(dirpath, "analysis_result.jld2")
    isfile(result_path) || error("$dirpath does not contain a analysis_result.jld2 file")

    return load(result_path, "analysis_result")
end

function summary(result::AbstractHLAAnalysisResult)
    df = DataFrame(name = String[], position = Int[], replacement = String[],
                   lower = Float64[], upper = Float64[], p_minus = Float64[])
    for model_result in result
        alleles = allele_posteriors(model_result)
        replacement = model_result.replacement

        for (allele, posterior) in alleles
            lower, upper = posterior_interval(posterior)
            p_minus = count(posterior .<= 0) / length(posterior)
            push!(df, [replacement.protein, replacement.position, 
                       string(replacement.replacement), lower, upper, p_minus])
        end
    end

    return df
end

function Base.iterate(result::AbstractHLAAnalysisResult, state = 1)
    files = result_files(result)
    state > length(files) && return nothing

    model_result = FileIO.load(files[state], "result")
    return (model_result, state + 1)
end

function result_files(result::AbstractHLAAnalysisResult)
    data_dirs = joinpath.(result.path, first(walkdir(result.path))[2])
    files = Vector{String}()

    for data_dir in data_dirs
        for file in readdir(data_dir)
            endswith(file, ".jld2") && push!(files, joinpath(data_dir, file))
        end
    end

    return files
end