@testset "binned_intervals(::AbstractVector{<: Real}, ::AbstractVector{<: Bool}" begin
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

    @test df[!, :expected][1] ≈ expected_1
    @test df[!, :expected][2] ≈ expected_2
    @test df[!, :observed][1] ≈ observed_1
    @test df[!, :observed][2] ≈ observed_2
    @test isapprox(df[!, :lower][1], interval_1[1], atol = 0.01)
    @test isapprox(df[!, :lower][2], interval_2[1], atol = 0.01)
    @test isapprox(df[!, :upper][1], interval_1[2], atol = 0.01)
    @test isapprox(df[!, :upper][2], interval_2[2], atol = 0.01)

    theta = range(0, 1, length = 100000)
    y = map(x -> rand(Bernoulli(x)), theta)
    df = Escape.binned_intervals(theta, y, bins = 10, lower = 0.0001, upper = 0.9999) 

    for row in eachrow(df) 
        @test row[:lower] <= row[:observed] <= row[:upper]
    end
end

@testset "check_calibration_arguments(::Any, ::Any)" begin
    @test Escape.check_calibration_arguments(
        [0.1, 0.2, 0.3], [true, false, true]
    ) === nothing
    @test_throws ErrorException Escape.check_calibration_arguments(
        ["a", "b", "c"], [true, false, true]
    )
    @test_throws ErrorException Escape.check_calibration_arguments(
        [0.1, 0.2], [1, 0]
    )
    @test_throws ErrorException Escape.check_calibration_arguments(
        [0.1, 0.2], [true, false, true]
    )
    @test_throws ErrorException Escape.check_calibration_arguments(
        [0.3, 1, 1.1], [true, false, true]
    )
end

@testset "@recipe function f(::Calibration)" begin
    theta = range(0, 1, length = 1000)
    y = map(x -> rand(Bernoulli(x)), theta)

    @test Escape.calibration(theta, y) isa Plots.Plot
end