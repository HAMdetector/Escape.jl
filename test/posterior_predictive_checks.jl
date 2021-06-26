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

@testset "@recipe f(::HLAModelResult)" begin
    result = @suppress Escape.run(
        Escape.HLAModel{1}(), Escape.HLADataset("Test").data[1], 
        iter = 10, chains = 1, warmup = 10
    )
    
    @test @suppress Escape.calibration_plot(result, samples = 10) isa Plots.Plot
end