using Test
using Unitful: Quantity, @u_str
using RNAstructure: efn2, energy

@testset "efn2" begin
    Tres = Tuple{Int, String, String, String}
    for (seq, dbns) in [
        "GGGAAACCC" => ["(((...)))", "((.....))"],
        ]
        for kwargs in [
            (; ),
            (; cmdline_opts="-s"),
            (; cmdline_opts=["-T", "300"]),
            (; cmdline_opts=["-T", 300]),
            ]
            dbn = first(dbns)
            @test efn2(seq, dbn; kwargs...) isa Tres
            @test efn2(seq, dbns; kwargs...) isa Tres
        end
    end
end

@testset "energy" begin
    Tres = typeof((0.0u"kcal/mol", 0.0u"kcal/mol"))
    for (seq, dbns) in [
        "GGGAAACCC" => ["(((...)))", "((.....))"],
        ]
        for kwargs in [
            (; ),
            (; cmdline_opts="-s"),
            (; cmdline_opts=["-T", "300"]),
            (; cmdline_opts=["-T", 300]),
            ]
            dbn = first(dbns)
            e = energy(seq, dbn; kwargs...)
            @test e isa Tres
            e = energy(seq, dbns; kwargs...)
            @test e isa Vector{Tres}
            @test length(e) == length(dbns)
        end
    end

    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        energy("", ""; cmdline_opts="-h")
    end

    # --writedetails option, parsing of detailed output not implemented
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        energy("", ""; cmdline_opts="-w")
    end
end
