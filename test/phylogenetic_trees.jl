@testset "extract_name(::String)" begin
    @test ismissing(Escape.extract_name(""))
    @test ismissing(Escape.extract_name("(A:3, B:1)"))
    @test ismissing(Escape.extract_name("(A:3, B:1):3"))
    @test ismissing(Escape.extract_name("(,)"))
    @test Escape.extract_name("B:5") == "B"
    @test Escape.extract_name("(A:3, B:1)E") == "E"
    @test Escape.extract_name("(A:3, B:1)E:5") == "E"
    @test Escape.extract_name("(,)E") == "E"
    @test Escape.extract_name("(,B)E:3") == "E"
end

@testset "extract_length(::String)" begin
    @test ismissing(Escape.extract_length(""))
    @test ismissing(Escape.extract_length("(A:3, B:1)"))
    @test ismissing(Escape.extract_length("(,)"))
    @test ismissing(Escape.extract_length("(A:3, B:1)E"))
    @test ismissing(Escape.extract_length("(,)E"))
    @test Escape.extract_length("(A:3, B:1):3") == 3.0
    @test Escape.extract_length("B:5") == 5.0
    @test Escape.extract_length("(A:3, B:1)E:5") == 5.0
    @test Escape.extract_length("(,B)E:3") == 3.0
end