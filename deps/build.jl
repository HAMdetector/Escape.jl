using StanInterface

for i in 1:4
    path = joinpath(@__DIR__, "..", "models", "model_$(i).stan")
    build_binary(path)
end