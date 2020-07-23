export HLAModelResultIO
export hla_model_result_io
export dir
export model
export dataset

struct HLAModelResultIO
    model::Escape.AbstractHLAModel
    ds::Escape.AbstractHLADataset
    dir::String
end

function Escape.run(
    model::Escape.AbstractHLAModel, 
    ds::Escape.AbstractHLADataset;
    result_dir::String = "",
    result_name::String = lowercase(split(string(model), '.')[2][1:end-2]),
    kwargs...
)   
    result_dir != "" || throw(ArgumentError("must specify keyword argument 'result_dir'"))
    isdir(result_dir) || throw(ArgumentError("directory $result_dir does not exist"))
    if result_name == "escape_result"
        throw(ArgumentError("$result_name must not be 'escape_result'"))
    end

    mkpath(joinpath(result_dir, result_name))
    for (i, data) in enumerate(ds.data)
        result_file = joinpath(result_dir, result_name, "$i.jls")

        if !isfile(result_file)
            res = Escape.run(model, data; kwargs...)
            serialize(joinpath(result_dir, result_name, "$i.jls"), res)
        end
    end

    result = HLAModelResultIO(model, ds, joinpath(result_dir, result_name))
    serialize(joinpath(result_dir, result_name, "result.jls"), result)

    return result
end

dir(result::HLAModelResultIO) = result.dir
model(result::HLAModelResultIO) = result.model
dataset(result::HLAModelResultIO) = result.ds

function hla_model_result_io(dir::String)
    result = deserialize(joinpath(dir, "result.jls"))

    return Escape.HLAModelResultIO(result.model, result.ds, dir)
end

function Base.getindex(result::HLAModelResultIO, i::Int)
    dir_ = dir(result)
    model_result = deserialize(joinpath(dir_, "$i.jls"))

    return model_result
end

function Base.iterate(iter::HLAModelResultIO)
    dir_ = dir(iter)
    model_result = deserialize(joinpath(dir_, "1.jls"))

    return (model_result, 1)
end

function Base.iterate(iter::HLAModelResultIO, state::Int)
    dir_ = dir(iter)
    model_path = joinpath(dir_, "$(state + 1).jls")

    if isfile(model_path)
        model_result = deserialize(joinpath(dir_, "$(state + 1).jls"))
        return (model_result, state + 1)
    else
        return nothing
    end
end

function Base.length(iter::HLAModelResultIO)
    return length(dataset(iter).data)
end