@testset "GGPlot(::RObject{VecSxp})" begin
    p = R"""
    suppressPackageStartupMessages(library("ggplot2"))
    df <- data.frame(a = rnorm(100), b = rnorm(100))
    p <- ggplot(df) +
        geom_point(aes(a, b))

    pdf(file = NULL)
    p
    """

    p = Escape.GGPlot(p)
    @test p isa Escape.GGPlot
end

@testset "save(::GGPlot)" begin
    p = R"""
    suppressPackageStartupMessages(library("ggplot2"))
    df <- data.frame(a = rnorm(100), b = rnorm(100))
    p <- ggplot(df) +
        geom_point(aes(a, b))

    pdf(file = NULL)
    p
    """

    p = Escape.GGPlot(p)
    temp_path = tempname() * "test.pdf"

    @test !isfile(temp_path)
    Escape.save(p, temp_path)
    @test isfile(temp_path)
    rm(temp_path)
end