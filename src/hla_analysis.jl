export HLAAnalysis, HLAAnalysisResult, analysis_result

@reexport using LightGraphs.SimpleGraphs
@reexport using MetaGraphs

struct HLAAnalysis{T <: HLAModel}
    name::String
    model::T
    dataset::AbstractHLADataset
end

struct HLAAnalysisResult{T <: HLAModel}
    analysis::HLAAnalysis{T}
    path::String
end

function run(analysis::HLAAnalysis{T}, dir::String; mincount::Int = 10) where T <: HLAModel
    isdir(dir) || error("$dir must be a directory.")
    to_filename(x::String) = replace(lowercase(x), ' ' => '_')

    root = joinpath(dir, to_filename(analysis.name))
    mkdir(root)
    analysis = rebase_analysis(analysis, root)

    for data in analysis.dataset.data
        data_dir = joinpath(root, to_filename(data.name))
        mkdir(data_dir)

        r = replacements(data, mincount = mincount)

        @sync @distributed for replacement in r
            result = run(analysis.model, data, replacement; wp = WorkerPool([myid()]))

            filename = string(replacement.protein, "_", replacement.position, 
                replacement.replacement, ".jld2")
            
            
            jldopen(joinpath(data_dir, filename), true, true, true, IOStream) do file
                file["result"] = result
            end
        end
    end

    analysis_result = HLAAnalysisResult{T}(analysis, root) 
    jldopen(joinpath(analysis_result.path, "analysis_result.jld2"), 
            true, true, true, IOStream) do file
        file["analysis_result"] = analysis_result
    end
end

function rebase_analysis(analysis::HLAAnalysis{T}, dir::String) where T <: HLAModel
    dataset = deepcopy(analysis.dataset)
    isdir(joinpath(dir, "data")) || mkdir(joinpath(dir, "data"))

    for data in analysis.dataset.data
        source = data.fasta_file
        destination = joinpath(dir, "data", basename(source))
        cp(source, destination, force = true)

        data.fasta_file = destination
    end

    return HLAAnalysis(analysis.name, analysis.model, dataset)
end

function analysis_result(dirpath::AbstractString)
    isdir(dirpath) || error("directory $dirpath does not exist.")
    result_path = joinpath(dirpath, "analysis_result.jld2")
    isfile(result_path) || error("$dirpath does not contain an analysis_result.jld2 file")
    
    result = JLD2.load(File(format"JLD2", result_path), "analysis_result")

    return HLAAnalysisResult(result.analysis, dirpath)
end

function summary(result::HLAAnalysisResult{T}) where T <: HLAModel
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

function summary(result::HLAAnalysisResult{FisherTest})
    df = DataFrame(name = String[], allele = String[], position = Int[],
                   replacement = String[], p = Float64[])
    
    p_values = Float64[]
    for (i, model_result) in enumerate(result)
        p_values = vcat(p_values, collect(values(model_result.p_values)))
    end
    sort!(p_values)

    p_adjusted = adjust(p_values, BenjaminiHochberg())
    threshold = p_values[findfirst(x -> x >= 0.05, p_adjusted)]

    for (i, model_result) in enumerate(result)
        alleles = keys(model_result.p_values)
        replacement = model_result.replacement

        for allele in alleles
            if model_result.p_values[allele] <= threshold
                push!(df, [replacement.protein, string(allele), replacement.position,
                           string(replacement.replacement), model_result.p_values[allele]])
            end
        end
    end

    return df
end

function Base.iterate(result::HLAAnalysisResult{T}, state = 1) where T <: HLAModel
    files = result_files(result)
    state > length(files) && return nothing

    model_result = FileIO.load(files[state], "result")
    return (model_result, state + 1)
end

function result_files(result::HLAAnalysisResult{T}) where T <: HLAModel
    data_dirs = joinpath.(result.path, first(walkdir(result.path))[2])
    files = Vector{String}()

    for data_dir in data_dirs
        for file in readdir(data_dir)
            endswith(file, ".jld2") && push!(files, joinpath(data_dir, file))
        end
    end

    return files
end