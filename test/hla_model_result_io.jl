let
    result_io = Escape.HLAModelResultIO(Escape.HLAModel{4}(), Escape.HLADataset("Test"), "a")

    @testset "HLAModelResultIO(::Escape.HLAModel, ::Escape.HLADataset)" begin
        @test result_io isa Escape.HLAModelResultIO
    end

    @testset "dir(::HLAModelResultIO)" begin
        @test Escape.dir(result_io) == "a"
    end

    @testset "model(::HLAModelResultIO)" begin
        @test Escape.model(result_io) == Escape.HLAModel{4}()
    end
 
    @testset "dataset(::HLAModelResultIO)" begin
        @test Escape.dataset(result_io) == result_io.ds
    end

    tmp = tempname()

    try
        result = @suppress Escape.run(
            Escape.HLAModel{4}(), Escape.HLADataset("Test"),
            mincount = 1,
            result_dir = tempdir(),
            result_name = splitdir(tmp)[end],
            iter = 10, warmup = 10, chains = 1
        )

        @testset "hla_model_result_io(::String)" begin
            @test Escape.hla_model_result_io(tmp) isa Escape.HLAModelResultIO
        end

        @testset "Base.getindex(::HLAModelResultIO, ::Int)" begin
            @test result[1] isa Escape.HLAModelResult
        end

        @testset "Base.iterate(::HLAModelResultIO)" begin
            i = 0
            for res in result
                i += 1
                @test res isa Escape.HLAModelResult
            end
            @test i == 2
        end

        @testset "Base.length(::HLAModelResultIO)" begin
            @test length(result) == 2
        end

    finally
        rm(tmp, force = true, recursive = true)
    end
end