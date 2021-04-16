@testset "phylogeny_information(::AbstractHLAModel, ::AbstractHLAData, ::Vector{Replacement})" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]
    r = Escape.replacements(ds.data[1])

    @test @suppress collect(Escape.phylogeny_information(m, data, r)) isa Matrix
end

@testset "epitope_information(::AbstractHLAModel, ::AbstractHLAData, ::Vector{Replacement}, ::Int)" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]
    r = Escape.replacements(ds.data[1])

    @test @suppress collect(Escape.epitope_information(m, data, r, 1)) isa Matrix
end

@testset "stan_input(::AbstractHLAModel, ::AbstractHLAData)" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]

    @test Escape.stan_input(m, data) isa Dict
end

@testset "Escape.run(::HLAModel, data::AbstractHLAData)" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]

    @test @suppress Escape.run(m, data, iter = 5, warmup = 5, chains = 1) isa 
        Escape.HLAModelResult
end

@testset "Escape.reduce_size!(::StanInterface.Stanfit)" begin
    ds = Escape.HLADataset("Test")
    m = Escape.HLAModel{4}()
    data = ds.data[1]

    res = @suppress Escape.run(m, data, iter = 5, warmup = 5, chains = 1, 
        keep_all_parameters = true)
    previous_size = Base.summarysize(res.sf)

    @test "aux1_tau.1" in keys(res.sf.result[1])
    Escape.reduce_size!(res.sf)

    new_size = Base.summarysize(res.sf)
    @test !("aux1_tau.1" in keys(res.sf.result[1]))
    @test new_size < previous_size
end