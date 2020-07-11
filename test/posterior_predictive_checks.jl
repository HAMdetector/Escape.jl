@testset "check_calibration_arguments(::AbstractVector{<: Real}, ::AbstractVector{<: Bool}))" begin
    @test isnothing(
        Escape.check_calibration_arguments([0.1, 0.2, 0.3], [true, false, true])
    )
    @test_throws ErrorException Escape.check_calibration_arguments(
        [0.1, 0.2], [true, false, true]
    )
    @test_throws ErrorException Escape.check_calibration_arguments(
        [0.3, 1, 1.1], [true, false, true]
    )
end

@testset "binned_intervals(::AbstractVector{<: Real}, ::AbstractVector{<: Bool})" begin
    theta = collect(0.3:0.1:1)
    y = [0, 0, 0, 1, 0, 1, 1, 1]

    expected_1 = mean(theta[1:4])
    expected_2 = mean(theta[5:8])
    observed_1 = mean(y[1:4])
    observed_2 = mean(y[5:8])

    interval_1 = quantile(rand(Beta(2, 4), 50000), (0.025, 0.975))
    interval_2 = quantile(rand(Beta(4, 2), 50000), (0.025, 0.975))

    @test Escape.binned_intervals(theta, Bool.(y), bins = 2) isa DataFrames.DataFrame
    df = Escape.binned_intervals(theta, Bool.(y), bins = 2)

    theta = range(0, 1, length = 100000)
    y = map(x -> rand(Bernoulli(x)), theta)
    df = Escape.binned_intervals(theta, y, bins = 10, lower = 0.0001, upper = 0.9999) 

    for row in eachrow(df) 
        @test row[:lower] <= row[:observed] <= row[:upper]
    end
end

@testset "@recipe f(::Type{Val{:calibration}}, x, y, z)" begin
    theta = range(0, 1, length = 1000)
    y = map(x -> rand(Bernoulli(x)), theta)

    @test Plots.plot(theta, y, seriestype = :calibration) isa Plots.Plot
    @test Plots.plot(theta, y, bins = 20, seriestype = :calibration) isa Plots.Plot
end

@testset "@recipe f(::HLAModelResult)" begin
    result = @suppress Escape.run(
        Escape.HLAModel{2}(), Escape.HLADataset("Test").data[1], 
        iter = 10, chains = 1, warmup = 10, mincount = 1
    )
    
    @test Escape.calibration_plot(result) isa Plots.Plot
end

@testset "@recipe f(::HLAModelResultIO)" begin
    tmp = tempname()

    try
        result = @suppress Escape.run(
            Escape.HLAModel{2}(), Escape.HLADataset("Test"),
            mincount = 1,
            result_dir = tempdir(),
            result_name = splitdir(tmp)[end],
            iter = 10, warmup = 10, chains = 1
        )

        @test Escape.calibration_plot(result) isa Plots.Plot
    finally
        rm(tmp, force = true, recursive = true)
    end
end