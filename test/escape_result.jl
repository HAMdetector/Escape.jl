@testset "dir(::EscapeResult)" begin
    result = Escape.EscapeResult(Escape.Model2(), Escape.HLADataset("Test"), "a")
    @test Escape.dir(result) == "a"
end

@testset "model(::EscapeResult)" begin
    result = Escape.EscapeResult(Escape.Model2(), Escape.HLADataset("Test"), "a")
    @test Escape.model(result) == Escape.Model2()
end
 
@testset "dataset(::EscapeResult)" begin
    result = Escape.EscapeResult(Escape.Model2(), Escape.HLADataset("Test"), "a")
    @test Escape.dataset(result) == result.ds
end

@testset "Base.iterate(::EscapeResult)" begin
    tmp = tempname()

    try
        result = @suppress Escape.run(
            Escape.Model2(), Escape.HLADataset("Test"), 
            mincount = 1,
            result_dir = tempdir(), 
            result_name = splitdir(tmp)[end],
            iter = 10, warmup = 10, chains = 1
        )

        i = 0
        for res in result
            i += 1
            @test res isa Escape.HLAModelResult
        end
        @test i == 2
    finally
        rm(tmp, force = true, recursive = true)
    end
end

@testset "Base.length(::EscapeResult)" begin
    result = Escape.EscapeResult(Escape.Model2(), Escape.HLADataset("Test"), "a")
    @test length(result) == 2
end