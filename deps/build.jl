using StanInterface

for i in 1:4
    path = joinpath(@__DIR__, "..", "models", "model_$(i).stan")
    isfile(joinpath(@__DIR__, "..", "models", "model_$i")) || build_binary(path)
end