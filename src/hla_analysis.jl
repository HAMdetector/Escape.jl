export HLAAnalysis, HLAAnalysisResult, analysis_result

struct HLAAnalysis{T <: HLAModel}
    name::String
    model::T
    dataset::AbstractHLADataset
end

struct HLAAnalysisResult{T <: HLAModel}
    analysis::HLAAnalysis{T}
    path::String
end

struct HLAAnalysisIterator
    files::Vector{String}
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
                replacement.replacement, ".jls")
            
            serialize(joinpath(data_dir, filename), result)
        end
    end

    analysis_result = HLAAnalysisResult{T}(analysis, root)
    serialize(joinpath(analysis_result.path, "analysis_result.jls"), analysis_result)

    return analysis_result
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
    result_path = joinpath(dirpath, "analysis_result.jls")
    isfile(result_path) || error("$dirpath does not contain an analysis_result.jls file")
    
    result = deserialize(result_path)

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

function summary(result::HLAAnalysisResult{FisherTest}; threshold::Float64 = 0.05)
    df = DataFrame(name = String[], allele = String[], position = Int[],
                   replacement = String[], log_odds = Float64[], p = Float64[],
                   p_adjusted = Float64[])
    
    p_values = Float64[]
    for (i, model_result) in enumerate(result)
        p_values = vcat(p_values, collect(values(model_result.p_values)))
    end
    sort!(p_values)

    p_adjusted = adjust(p_values, BenjaminiHochberg())
    d = Dict(zip(p_values, p_adjusted))

    for (i, model_result) in enumerate(result)
        alleles = keys(model_result.p_values)
        replacement = model_result.replacement

        for allele in alleles
            if model_result.p_values[allele] <= 0.05
                push!(df, [replacement.protein, string(allele), replacement.position,
                           string(replacement.replacement), model_result.log_odds[allele],
                           model_result.p_values[allele], d[model_result.p_values[allele]]])
            end
        end
    end

    return df
end


function result_files(result::HLAAnalysisResult{T}) where T <: HLAModel
    data_dirs = joinpath.(result.path, first(walkdir(result.path))[2])
    files = Vector{String}()

    for data_dir in data_dirs
        for file in readdir(data_dir)
            endswith(file, ".jls") && push!(files, joinpath(data_dir, file))
        end
    end

    return files
end

function Base.iterate(I::HLAAnalysisIterator, state = 1) where T <: HLAModel
    state > length(I.files) && return nothing

    model_result = deserialize(I.files[state])
    return (model_result, state + 1)
end

function Base.getindex(I::HLAAnalysisIterator, i::Int)
    return deserialize(I.files[i])
end

function Base.getindex(result::HLAAnalysisResult, collection)
    files = result_files(result)

    return HLAAnalysisIterator(files[collection])
end

function Base.iterate(result::HLAAnalysisResult, state = 1)
    files = result_files(result)
    state > length(files) && return nothing
    
    return (HLAAnalysisIterator(files)[state], state + 1)
end