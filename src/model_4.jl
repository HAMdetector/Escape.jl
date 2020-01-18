function stan_input(
    model::Model4, data::AbstractHLAData; 
    depth::Int = 1, mincount::Int = 10
)
    X = hla_matrix(data.hla_types; depth = depth)
    r = replacements(data, mincount = mincount)
    N = size(X)[1]
    D = size(X)[2]
    R = length(r)

    ys = Matrix{Int64}(undef, R, N + 1)
    xs = Matrix{Float64}(undef, R, N * (D + 1))
    y_mean = Vector{Float64}(undef, R)
    hla_mean = [mean(X[:, i]) for i in 1:D]
    fill!(ys, -1)
    fill!(xs, -1)

    phy = @distributed hcat for i in 1:R
        phylogeny_probabilities(r[i], data, data.tree)
    end

    for i in 1:R
        t = targets(r[i], data)
        y_mean[i] = mean(skipmissing(t))
        count = 0
        for j in 1:N
            ys[i, ]
            if !ismissing(t[j])
                count += 1
                ys[i, count + 1] = t[j]
            end
        end

        ys[i, 1] = count
        X_ = hcat(phy[.!ismissing.(t), i], X[.!ismissing.(t), :])
        xs[i, 1:((D + 1) * count)] = reshape(X_', :, 1)
    end


    d = Dict("N" => N, "D" => D, "R" => R, "ys" => ys, "xs" => xs,
             "y_mean" => y_mean, "hla_mean" => hla_mean, "p0" => 5)
    
    return d
end