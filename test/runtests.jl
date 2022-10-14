using Test
using Unitful: Quantity, @u_str
using RNAstructure: energy

@testset "energy" begin
    Tres = typeof((0.0u"kcal/mol", 0.0u"kcal/mol"))
    @test energy("GGGAAACCC", "(((...)))") isa Tres
    @test energy("GGGAAACCC", "(((...)))"; cmdline_opts=["-T", "300"]) isa Tres

    Tres = Vector{Tres}
    e = energy("GGGAAACCC", ["(((...)))"])
    @test e isa Tres
    @test length(e) == 1
    e = energy("GGGAAACCC", ["(((...)))"]; cmdline_opts=["-T", "300"])
    @test e isa Tres
    @test length(e) == 1
    e = energy("GGGAAACCC", ["(((...)))", "((.....))"])
    @test e isa Tres
    @test length(e) == 2
    e = energy("GGGAAACCC", ["(((...)))", "((.....))"]; cmdline_opts=["-T", "300"])
    @test e isa Tres
    @test length(e) == 2
end
