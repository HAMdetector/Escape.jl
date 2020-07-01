@testset "HLAModelResultIO(::Escape.HLAModel, ::Escape.HLADataset)" begin
    ds = Escape.HLADataset("Rousseau")
    @test Escape.HLAModelResultIO(Escape.HLAModel{4}(), ds, "") isa Escape.HLAModelResultIO
end

@testset "dir(::HLAModelResultIO)" begin
    result = Escape.HLAModelResultIO(Escape.HLAModel{4}(), Escape.HLADataset("Test"), "a")
    @test Escape.dir(result) == "a"
end

@testset "model(::HLAModelResultIO)" begin
    result = Escape.HLAModelResultIO(Escape.HLAModel{4}(), Escape.HLADataset("Test"), "a")
    @test Escape.model(result) == Escape.HLAModel{4}()
end
 
@testset "dataset(::HLAModelResultIO)" begin
    result = Escape.HLAModelResultIO(Escape.HLAModel{4}(), Escape.HLADataset("Test"), "a")
    @test Escape.dataset(result) == result.ds
end

@testset "hla_model_result_io(::String)" begin
    tmp = tempname()

    try
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test"),
            mincount = 1,
            result_dir = tempdir(),
            result_name = splitdir(tmp)[end],
            iter = 10, warmup = 10, chains = 1
        )

        @test Escape.hla_model_result_io(tmp) isa Escape.HLAModelResultIO
    finally
        rm(tmp, force = true, recursive = true)
    end
end

@testset "Base.getindex(::HLAModelResultIO, ::Int)" begin
    tmp = tempname()

    try
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test"),
            mincount = 1,
            result_dir = tempdir(),
            result_name = splitdir(tmp)[end],
            iter = 10, warmup = 10, chains = 1
        )

        @test result[1] isa Escape.HLAModelResult
    finally
        rm(tmp, force = true, recursive = true)
    end
end

@testset "Base.iterate(::HLAModelResultIO)" begin
    tmp = tempname()

    try
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test"), 
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

@testset "Base.length(::HLAModelResultIO)" begin
    result = Escape.HLAModelResultIO(Escape.HLAModel{4}(), Escape.HLADataset("Test"), "a")
    @test length(result) == 2
end