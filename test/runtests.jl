using Test
using Unitful: Quantity, @u_str
using RNAstructure: efn2, energy, design, fold, edcalculator, ensemble_defect

@testset "efn2" begin
    Tres = Tuple{Int, String, String, String}
    for (seq, dbns) in [
        "GGGAAACCC" => ["(((...)))", "((.....))"],
        ]
        for kwargs in [
            (; ),
            (; cmdargs="-s"),
            (; cmdargs=`-T 300`),
            (; cmdargs=["-T", 300]),
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
            (; cmdargs="-s"),
            (; cmdargs=`-T 300`),
            (; cmdargs=["-T", 300]),
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
        energy("", ""; cmdargs=`-h`)
    end

    # --writedetails option, parsing of detailed output not implemented
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        energy("", ""; cmdargs=`-w`)
    end
end

@testset "design" begin
    Tres = typeof((; seq = "", seed = ""))
    target = "(((...)))"
    for kwargs in [
        (; ),
        (; cmdargs=`-s 42`),
        ]
        res = design(target; kwargs...)
        @test res isa Tres
        @test length(res.seq) == length(target)
    end

    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        design(""; cmdargs=`-h`)
    end
end

@testset "fold" begin
    Tres = Tuple{Int,String,String,String}
    seq = "GGGAAAACCC"

    for kwargs in [
        (; ),
        (; cmdargs=`-h`),
        (; cmdargs=`-mfe`),
        ]
        res = fold(seq; kwargs...)
        @test res isa Tres
    end
end

@testset "edcalculator" begin
    Tres = Tuple{Int,String,String}
    seq = "GGGAAAACCC"
    dbn = "(((....)))"
    for kwargs in [
        (; ),
        (; cmdargs=`-h`),
        (; cmdargs=`-s 1`),
        ]
        res = edcalculator(seq, dbn; kwargs...)
        @test res isa Tres
        res = edcalculator(seq, [dbn, dbn]; kwargs...)
        @test res isa Tres
    end
end

@testset "ensemble_defect" begin
    Tres = Tuple{Float64,Float64}
    seq = "GGGAAAACCC"
    dbn = "(((....)))"
    for kwargs in [
        (; ),
        (; cmdargs=`-s 1`),
        ]
        res = ensemble_defect(seq, dbn; kwargs...)
        @test res isa Tres
        res = ensemble_defect(seq, [dbn]; kwargs...)
        @test res isa Vector{Tres}
        res = ensemble_defect(seq, [dbn, dbn]; kwargs...)
        @test res isa Vector{Tres}
    end

    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        ensemble_defect("", ""; cmdargs=`-h`)
    end

end
