using RNAstructure: Plot

@testset "plot" begin
    Tres = Plot
    dbn = "(((...)))"
    seq = "GGGAAACCC"
    for kwargs in [
        (; ),
        (; args=`--flat`),
        (; args=`--circle`),
        (; args=`--uncircled`),
        (; args=`--levorotatory`),
        ]
        p = plot(dbn, seq; kwargs...)
        @test p isa Tres
        @test length(p.svg) > 0
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        plot("", ""; args=`-h`)
    end

end
