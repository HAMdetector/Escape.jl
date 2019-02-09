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

function run(analysis::HLAAnalysis, dir::String)
    isdir(dir) || error("$dir must be a directory.")
    to_filename(x::String) = replace(lowercase(x), ' ' => '_')

    root = joinpath(dir, to_filename(analysis.name))
    mkdir(root)
    analysis = rebase_analysis(analysis, root)

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
            
            
            jldopen(joinpath(data_dir, filename), true, true, true, IOStream) do file
                file["result"] = result
            end
        end
    end

    analysis_result = HLAAnalysisResult(analysis, root) 
    jldopen(joinpath(analysis_result.path, "analysis_result.jld2"), 
            true, true, true, IOStream) do file
        file["analysis_result"] = analysis_result
    end
end

function rebase_analysis(analysis::HLAAnalysis, dir::String)
    dataset = deepcopy(analysis.dataset)
    isdir(joinpath(dir, "data")) || mkdir(joinpath(dir, "data"))

    for data in analysis.dataset.data
        source = data.fasta_file.path
        destination = joinpath(dir, "data", basename(source))
        cp(source, destination, force = true)

        data.fasta_file.path = destination
    end

    return HLAAnalysis(analysis.name, analysis.model, dataset)
end

function analysis_result(dirpath::AbstractString)
    isdir(dirpath) || error("directory $dirpath does not exist.")
    result_path = joinpath(dirpath, "analysis_result.jld2")
    isfile(result_path) || error("$dirpath does not contain an analysis_result.jld2 file")

    return load(result_path, "analysis_result")
end

function summary(result::AbstractHLAAnalysisResult)
    df = DataFrame(name = String[], allele = String[], position = Int[], 
                   replacement = String[], lower = Float64[], upper = Float64[], 
                   p_minus = Float64[])

    for (i, model_result) in enumerate(result)
        alleles = relevant_alleles(model_result)
        replacement = model_result.replacement

        for (allele, posterior) in alleles
            lower, upper = posterior_interval(posterior)
            p_minus = count(posterior .<= 0) / length(posterior)
            push!(df, [replacement.protein, string(allele), replacement.position, 
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