using Test
using Unitful: Quantity, @u_str
using RNAstructure
using RNAstructure: run_draw, run_EDcalculator, run_efn2,
    run_EnsembleEnergy, run_Fold, run_MaxExpect, run_partition!,
    run_ProbabilityPlot, run_stochastic

include("parse-ct-format.jl")

@testset "bpp" begin
    Tres = Matrix{Float64}
    for seq in ["GGGAAACCC", "AAAAAAA"]
        for kwargs in [
            (; ),
            (; cmdargs=`-T 300`),
            ]
            n = length(seq)
            p = bpp(seq; kwargs...)
            @test p isa Tres
            @test size(p) == (n,n)
            @test all(x -> 0 <= x <= 1, p)
        end
    end

    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        bpp(""; cmdargs=`-h`)
    end
end

@testset "design" begin
    Tres = typeof((; seq = "", seed = ""))
    target = "(((...)))"
    for kwargs in [
        (; ),
        (; cmdargs=`--DNA`),
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

    # --DNA doesn't work yet (because for DNA no uncertainties are returned)
    @test_broken redirect_stdio(stdout=devnull, stderr=devnull) do
        energy("GGGAAACCC", "(((...)))"; cmdargs=`--DNA`)
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

@testset "mea" begin
    Tres = Vector{String}
    for seq in ["GGGAAACCC", "AAAAAAA"]
        for kwargs in [
            (; ),
            (; cmdargs_partition=`-T 300`),
            (; cmdargs_maxexpect=`-w 0`),
            ]
            n = length(seq)
            res = mea(seq; kwargs...)
            @test res isa Tres
            @test all(dbn -> length(dbn) == n, res)
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        mea(""; cmdargs_partition=`-h`)
    end
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        mea(""; cmdargs_maxexpect=`-h`)
    end
end

@testset "mfe" begin
    Tres = Tuple{typeof(0.0u"kcal/mol"),String}
    for seq in ["GGGAAAACCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; cmdargs=`-T 300`),
            ]
            res = mfe(seq; kwargs...)
            @test res isa Tres
        end
    end
end

@testset "partfn" begin
    Tres = typeof(0.0u"kcal/mol")
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; cmdargs=`--DNA`),
            ]
            res = partfn(seq; kwargs...)
            @test res isa Tres
        end
    end
end

@testset "prob_of_structure" begin
    Tres = Float64
    for inputs in [
        ("GGGAAACCC",
         "(((...)))"),
        ("AAAAAAAAA",
         ".........")]
        for kwargs in [
            (; ),
            (; cmdargs=``),  # TODO: no common cmdargs yet
            ]
            res = prob_of_structure(inputs...; kwargs...)
            @test res isa Tres
        end
    end

    # --DNA doesn't work yet (because for DNA no uncertainties are returned)
    @test_broken redirect_stdio(stdout=devnull, stderr=devnull) do
        prob_of_structure("GGGAAACCC", "(((...)))"; cmdargs=`--DNA`)
    end
end

@testset "sample_structures" begin
    Tres = Vector{String}
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; cmdargs=`-e 100`),
            ]
            res = sample_structures(seq; kwargs...)
            @test res isa Tres
            @test length(res) > 0
            @test all(dbn -> length(dbn) == length(seq), res)
        end
    end
end

@testset "subopt" begin
    Tres = Vector{Tuple{String,typeof(0.0u"kcal/mol")}}
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; cmdargs=`-T 300`),
            ]
            res = subopt(seq; kwargs...)
            @test res isa Tres
            @test length(res) > 0
            @test all(dbn_en -> length(dbn_en[1]) == length(seq), res)
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        subopt(""; cmdargs=`-h`)
    end
end

@testset "run_draw" begin
    Tres = Tuple{Int,String,String,String}
    inputdata = [
        ("(((...)))",
         "GGGAAACCC"),
        (".........",
         "NNNNNNNNN"),
    ]
    for inputs in inputdata
        for kwargs in [
            (; ),
            (; cmdargs=`-l`),
            (; cmdargs=`--svg`),
            ]
            res = run_draw(inputs...; kwargs...)
            @test res isa Tres
        end
    end
end


@testset "run_EDcalculator" begin
    Tres = Tuple{Int,String,String}
    seq = "GGGAAAACCC"
    dbn = "(((....)))"
    for kwargs in [
        (; ),
        (; cmdargs=`-h`),
        (; cmdargs=`-s 1`),
        ]
        res = run_EDcalculator(seq, dbn; kwargs...)
        @test res isa Tres
        res = run_EDcalculator(seq, [dbn, dbn]; kwargs...)
        @test res isa Tres
    end
end

@testset "run_efn2" begin
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
            @test run_efn2(seq, dbn; kwargs...) isa Tres
            @test run_efn2(seq, dbns; kwargs...) isa Tres
        end
    end
end

@testset "run_EnsembleEnergy" begin
    Tres = Tuple{Int,String,String}
    for seq in [
        "GGGAAACCC",
        ]
        for kwargs in [
            (; ),
            (; cmdargs="--DNA"),
            (; cmdargs=`-T 300`),
            ]
            @test run_EnsembleEnergy(seq; kwargs...) isa Tres
        end
    end
end

@testset "run_Fold" begin
    Tres = Tuple{Int,String,String,String}
    seq = "GGGAAAACCC"

    for kwargs in [
        (; ),
        (; cmdargs=`-h`),
        (; cmdargs=`-mfe`),
        ]
        res = run_Fold(seq; kwargs...)
        @test res isa Tres
    end
end

@testset "run_MaxExpect" begin
    Tres = Tuple{Int,String,String,String}
    for seq in [
        "GGGAAAACCC",
        "AAAAAAA",
        ]
        for kwargs in [
            (; ),
            (; cmdargs=`-p 25`),
            ]
            mktemp() do pf_savefile, _
                ps, out, err = run_partition!(pf_savefile, seq)
                if ps != 0
                    @warn "run_partition returned non-zero exit status"
                end
                res = run_MaxExpect(pf_savefile; kwargs...)
                @test res isa Tres
            end
        end
    end
end

@testset "run_partition!" begin
    Tres = Tuple{Int,String,String}
    seq = "GGGAAAACCC"

    for kwargs in [
        (; ),
        (; cmdargs=`-T 300`),
        ]
        mktemp() do pf_savefile, _
            res = run_partition!(pf_savefile, seq; kwargs...)
            @test res isa Tres
        end
    end
end

@testset "run_ProbabilityPlot" begin
    Tres = Tuple{Int,String,String,String}
    seq = "GGGAAAACCC"

    for kwargs in [
        (; ),
        (; cmdargs=`-min 0.1`),
        ]
        mktemp() do pf_savefile, _
            ps, out, err = run_partition!(pf_savefile, seq)
            if ps != 0
                @warn "run_partition returned non-zero exit status"
            end
            res = run_ProbabilityPlot(pf_savefile; kwargs...)
            @test res isa Tres
        end
    end
end

@testset "run_stochastic" begin
    Tres = Tuple{Int,String,String,String}
    seq = "GGGAAAACCC"

    for kwargs in [
        (; ),
        (; cmdargs=`-s 42`),
        ]
        res = run_stochastic(seq; kwargs...)
        @test res isa Tres
    end
end
