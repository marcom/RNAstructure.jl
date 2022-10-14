using Test
using Unitful: Quantity
using RNAstructure: energy

@testset "energy" begin
    @test energy("GGGAAACCC", "(((...)))") isa Tuple{Quantity, Quantity}
    @test energy("GGGAAACCC", "(((...)))"; cmdline_opts=["-T", "300"]) isa Tuple{Quantity, Quantity}
end
