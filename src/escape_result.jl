export EscapeResult
export escape_result
export dir
export model
export dataset

struct EscapeResult
    model::Escape.HLAModel
    ds::Escape.AbstractHLADataset
    dir::String
end

function Escape.run(
    model::Escape.HLAModel, 
    ds::Escape.AbstractHLADataset;
    result_dir::String = "",
    result_name::String = lowercase(split(string(model), '.')[2][1:end-2]),
    kwargs...
)   
    result_dir != "" || throw(ArgumentError("must specify keyword argument 'result_dir'"))
    isdir(result_dir) || throw(ArgumentError("directory $result_dir does not exist"))
    if result_name == "escape_result"
        throw(ArgumentError("$result_name must not be 'escape_result'"))
    elseif isdir(joinpath(result_dir, result_name))
        throw(ArgumentError("a directory with name $result_name already exists."))
    end

    mkpath(joinpath(result_dir, result_name))
    for (i, data) in enumerate(ds.data)
        res = Escape.run(model, data; kwargs...)

        serialize(joinpath(result_dir, result_name, "$i.jls"), res)
    end

    result = EscapeResult(model, ds, joinpath(result_dir, result_name))
    serialize(joinpath(result_dir, result_name, "result.jls"), result)

    return result
end

dir(result::EscapeResult) = result.dir
model(result::EscapeResult) = result.model
dataset(result::EscapeResult) = result.ds

function escape_result(dir::String)
    result = deserialize(joinpath(dir, "result.jls"))

    return result
end

function Base.getindex(result::EscapeResult, i::Int)
    dir_ = dir(result)
    model_result = deserialize(joinpath(dir_, "$i.jls"))

    return model_result
end

function Base.iterate(iter::EscapeResult)
    dir_ = dir(iter)
    model_result = deserialize(joinpath(dir_, "1.jls"))

    return (model_result, 1)
end

function Base.iterate(iter::EscapeResult, state::Int)
    dir_ = dir(iter)
    model_path = joinpath(dir_, "$(state + 1).jls")

    if isfile(model_path)
        model_result = deserialize(joinpath(dir_, "$(state + 1).jls"))
        return (model_result, state + 1)
    else
        return nothing
    end
end

function Base.length(iter::EscapeResult)
    return length(dataset(iter).data)
end