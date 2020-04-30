using StanInterface

for i in 1:4
    path = joinpath(@__DIR__, "..", "data", "stan", "model_$(i).stan")
    build_binary(path)
end