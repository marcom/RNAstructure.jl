using Test
using Unitful: Quantity, @u_str
using Random: randstring
using RNAstructure
using RNAstructure: run_AllSub, run_ct2dot, run_CycleFold, run_draw,
    run_dot2ct, run_EDcalculator, run_efn2, run_EnsembleEnergy,
    run_Fold, run_MaxExpect, run_partition!, run_ProbabilityPlot,
    run_RemovePseudoknots, run_stochastic

include("ct-format.jl")
include("plot.jl")

const DBN_CT = [
    (; title = "",
     seq = "NNNNNNNNN",
     dbns = ["(((...)))"]) =>
         """
         9
         1 N       0    2    9    1
         2 N       1    3    8    2
         3 N       2    4    7    3
         4 N       3    5    0    4
         5 N       4    6    0    5
         6 N       5    7    0    6
         7 N       6    8    3    7
         8 N       7    9    2    8
         9 N       8    0    1    9
         """,
    (; title = "Fooo bar",
     seq = "GGGGAAUCCCC",
     dbns = ["((.(...)))."]) =>
         """
         11 Fooo bar
          1 G       0    2   10    1
          2 G       1    3    9    2
          3 G       2    4    0    3
          4 G       3    5    8    4
          5 A       4    6    0    5
          6 A       5    7    0    6
          7 U       6    8    0    7
          8 C       7    9    4    8
          9 C       8   10    2    9
         10 C       9   11    1   10
         11 C      10    0    0   11
         """,
]

@testset "bpp" begin
    Tres = Matrix{Float64}
    for seq in ["GGGAAACCC", "AAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`-T 300`),
            ]
            n = length(seq)
            p = bpp(seq; kwargs...)
            @test p isa Tres
            @test size(p) == (n,n)
            @test all(x -> 0 <= x <= 1, p)
            @test all(I -> ((i,j) = Tuple(I); p[i,j] == p[j,i]), CartesianIndices(p))
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        bpp(""; args=`-h`)
    end
end

@testset "ct2dbn" begin
    Tres = typeof((; title="", seq="", dbn=""))

    ct = first(DBN_CT)[2]
    @test ct2dbn(ct) isa Tres

    for (edbn, ct) in DBN_CT
        for i = 1:length(edbn.dbns)
            res = ct2dbn(ct, i)
            @test res isa Tres
            (; title, seq, dbn) = res
            @test seq == edbn.seq
            @test dbn == edbn.dbns[i]
        end
    end

    # broken: converting one of multiple structures if they have
    # different length
    # Known broken: RNAstructure-6.4.0
    ct1 = dbn2ct("(((...)))"; title="A structure", seq="NNNAAANNN")
    ct2 = dbn2ct("(...)"; title="A shorter structure", seq="AANNN")
    ct = ct1 * ct2
    @test_broken ct2dbn(ct)
end

@testset "cyclefold_bpp" begin
    Tres = Matrix{Float64}
    for seq in ["GGGAAAACCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`--bigloops`),
            ]
            n = length(seq)
            pij = cyclefold_bpp(seq; kwargs...)
            @test pij isa Tres
            @test size(pij) == (n, n)
            @test all(x -> 0.0 <= x <= 1.0, pij)
        end
    end
    # -h, --help option
    @test_throws ErrorException cyclefold_bpp(""; args=`-h`)
    @test_throws ErrorException cyclefold_bpp(""; args=`--help`)
end

@testset "cyclefold_mea" begin
    Tres = Vector{Int}
    for seq in ["GGGAAAACCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`--bigloops`),
            ]
            pt = cyclefold_mea(seq; kwargs...)
            @test pt isa Tres
            @test length(pt) == length(seq)
            @test all(i -> i >= 0, pt)
        end
    end
    # -h, --help option
    @test_throws ErrorException cyclefold_mea(""; args=`-h`)
    @test_throws ErrorException cyclefold_mea(""; args=`--help`)
end

@testset "cyclefold_mfe" begin
    Tres = Tuple{typeof(0.0u"kJ/mol"),Vector{Int}}
    for seq in ["GGGAAAACCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`--bigloops`),
            ]
            res = cyclefold_mfe(seq; kwargs...)
            @test res isa Tres
            en, pt = res
            @test length(pt) == length(seq)
            @test all(i -> i >= 0, pt)
        end
    end
    # -h, --help option
    @test_throws ErrorException cyclefold_mfe(""; args=`-h`)
    @test_throws ErrorException cyclefold_mfe(""; args=`--help`)
end

@testset "dbn2ct" begin
    Tres = String
    input_dbns = [
        "(((...)))",
        "......",
        "(((...[[[...)))...]]]",
        "(((...[[[...{{{...)))...]]]...}}}",
        ["(...)", "....."],
    ]
    for dbn in input_dbns
        res = dbn2ct(dbn)
        @test res isa Tres
        @test length(res) > 0

        title = "Foo bar"
        n = if dbn isa Vector
            length(first(dbn))
        else
            length(dbn)
        end
        seq = randstring("ACGU", n)
        res = dbn2ct(dbn; title, seq)
        @test res isa Tres
        @test length(res) > 0
    end

    # throw when structures have different length
    @test_throws ArgumentError dbn2ct([".", ".."])
    # throw when sequence has different length than structures
    @test_throws ArgumentError dbn2ct(".."; seq="N")
    @test_throws ArgumentError dbn2ct(["()", ".."]; seq="N")
end

@testset "dbn2ct |> ct2dbn" begin
    title_seq_dbn = [
        (; title = "Foo bar", seq = "GGUAAAACC", dbn = "(((...)))"),
    ]
    for (; title, seq, dbn) in title_seq_dbn
        edbn = dbn2ct(dbn; title, seq) |> ct2dbn
        @test (edbn.title, edbn.seq, edbn.dbn) == (title, seq, dbn)
    end

    # multiple structures
    input = (; title="Foo bar", seq="NN", dbns=["()", ".."])
    ct = dbn2ct(input.dbns; title=input.title, seq=input.seq)
    for i = 1:length(input.dbns)
        (; title, seq, dbn) = ct2dbn(ct, i)
        @test title == input.title
        @test seq == input.seq
        @test dbn == input.dbns[i]
    end

    # mutiple structures roundtrip doesn't work
    @test_broken dbn2ct(["()", ".."]; title="Foo bar", seq="NN") |> ct2dbn ==
        (; title="Foo bar", seq="NN", dbns=["()", ".."])
end

@testset "design" begin
    Tres = typeof((; seq = "", seed = ""))
    target = "(((...)))"
    for kwargs in [
        (; ),
        (; args=`-s 42`),
        ]
        res = design(target; kwargs...)
        @test res isa Tres
        @test length(res.seq) == length(target)
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        design(""; args=`-h`)
    end
end

@testset "energy" begin
    Tres = typeof((0.0u"kcal/mol", 0.0u"kcal/mol"))
    for (seq, dbns) in [
        "GGGAAACCC" => ["(((...)))", "((.....))"],
        ]
        for kwargs in [
            (; ),
            (; args="-s"),
            (; args=`-T 300`),
            (; args=["-T", 300]),
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
        energy("GGGAAACCC", "(((...)))"; args=`--DNA`)
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        energy("", ""; args=`-h`)
    end
    # --writedetails option, parsing of detailed output not implemented
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        energy("", ""; args=`-w`)
    end
end

@testset "ensemble_defect" begin
    Tres = Tuple{Float64,Float64}
    seq = "GGGAAAACCC"
    dbn = "(((....)))"
    for kwargs in [
        (; ),
        (; args=`-s 1`),
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
        ensemble_defect("", ""; args=`-h`)
    end
end

@testset "mea" begin
    Tres = Vector{String}
    for seq in ["GGGAAACCC", "AAAAAAA"]
        for kwargs in [
            (; ),
            (; args_partition=`-T 300`),
            (; args_maxexpect=`-w 0`),
            ]
            n = length(seq)
            res = mea(seq; kwargs...)
            @test res isa Tres
            @test all(dbn -> length(dbn) == n, res)
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        mea(""; args_partition=`-h`)
    end
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        mea(""; args_maxexpect=`-h`)
    end
end

@testset "mfe" begin
    Tres = Tuple{typeof(0.0u"kcal/mol"),String}
    for seq in ["GGGAAAACCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`-T 300`),
            ]
            res = mfe(seq; kwargs...)
            @test res isa Tres
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        mfe(""; args=`-h`)
    end
end

@testset "partfn" begin
    Tres = typeof(0.0u"kcal/mol")
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`--DNA`),
            ]
            res = partfn(seq; kwargs...)
            @test res isa Tres
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        partfn(""; args=`-h`)
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
            (; args=``),  # TODO: no common args yet
            ]
            res = prob_of_structure(inputs...; kwargs...)
            @test res isa Tres
        end
    end
    # --DNA doesn't work yet (because for DNA no uncertainties are returned)
    @test_broken redirect_stdio(stdout=devnull, stderr=devnull) do
        prob_of_structure("GGGAAACCC", "(((...)))"; args=`--DNA`)
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        prob_of_structure("", ""; args=`-h`)
    end
end

@testset "remove_pknots" begin
    dbns = [
        "......",
        "(((...)))",
        "(((...[[[...)))...]]]",
    ]
    for dbn in dbns
        pkfree_dbn = remove_pknots(dbn)
        @test length(pkfree_dbn) == length(dbn)
        # TODO: can we test more?
        # - pknot-free
        # - some basepairs of original structure are still there
    end
end

@testset "sample_structures" begin
    Tres = Vector{String}
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`-e 100`),
            ]
            res = sample_structures(seq; kwargs...)
            @test res isa Tres
            @test length(res) > 0
            @test all(dbn -> length(dbn) == length(seq), res)
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        sample_structures(""; args=`-h`)
    end
end

@testset "subopt" begin
    Tres = Vector{Tuple{String,typeof(0.0u"kcal/mol")}}
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`-T 300`),
            ]
            res = subopt(seq; kwargs...)
            @test res isa Tres
            @test length(res) > 0
            @test all(dbn_en -> length(dbn_en[1]) == length(seq), res)
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        subopt(""; args=`-h`)
    end
end

@testset "subopt_all" begin
    Tres = Vector{Tuple{String,typeof(0.0u"kcal/mol")}}
    for seq in ["GGGGAAACCCC", "AAAAAAAAAA"]
        for kwargs in [
            (; ),
            (; args=`-T 300`),
            ]
            res = subopt_all(seq; kwargs...)
            @test res isa Tres
            @test length(res) > 0
            @test all(dbn_en -> length(dbn_en[1]) == length(seq), res)
        end
    end
    # --help option
    @test_throws ErrorException redirect_stdio(stdout=devnull, stderr=devnull) do
        subopt_all(""; args=`-h`)
    end
end

@testset "run_AllSub" begin
    Tres = Tuple{Int,String,String,String}
    inputs = [
        "GGGAAACCC",
        "AAAAAAA",
    ]
    for seq in inputs
        for kwargs in [
            (; ),
            (; args=`-T 300`),
            (; args=`-p 100 -a 1000`),
            ]
            res = run_AllSub(seq; kwargs...)
            @test res isa Tres
        end
    end
end

@testset "run_ct2dot" begin
    Tres = Tuple{Int,String,String,String}
    for (_, ct) in DBN_CT
        for kwargs in [
            (; ),
            (; args=`-h`),
            ]
            res = run_ct2dot(ct; kwargs...)
            @test res isa Tres
        end
    end
end

@testset "run_CycleFold" begin
    Tres = Tuple{Int,String,String}
    inputs = [
        "GGGAAACCC",
        "AAAAAAA",
        ["GGGAAACCC", "AAAAAAA"],
    ]
    for input in inputs
        for kwargs in [
            (; ),
            (; args=`-p`),
            (; args=`-m`),
            (; args=`--turbo`),
            (; args=`--turbo -p`),
            ]
            res = run_CycleFold(input; kwargs...)
            @test res isa Tres
            ps, out, err = res
            @test ps == 0
            @test length(err) == 0
            @test length(out) > 0
        end
    end
end

@testset "run_dot2ct" begin
    Tres = Tuple{Int,String,String}
    input_dbns = [
        "(((...)))",
        "(((...[[[...)))...]]]",
        "(((...)))",
        ".........",
        ["(...)", "....."],
    ]
    for dbn in input_dbns
        for kwargs in [
            (; ),
            (; args=``),
            ]
            n = if dbn isa Vector
                length(first(dbn))
            else
                length(dbn)
            end
            seq = randstring("ACGU", n)

            res = run_dot2ct(dbn; kwargs...)
            @test res isa Tres
            res = run_dot2ct(dbn; seq, kwargs...)
            @test res isa Tres
        end
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
            (; args=`-l`),
            (; args=`--svg`),
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
        (; args=`-h`),
        (; args=`-s 1`),
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
            (; args="-s"),
            (; args=`-T 300`),
            (; args=["-T", 300]),
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
            (; args="--DNA"),
            (; args=`-T 300`),
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
        (; args=`-h`),
        (; args=`-mfe`),
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
            (; args=`-p 25`),
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
        (; args=`-T 300`),
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
        (; args=`-min 0.1`),
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

@testset "run_RemovePseudoknots" begin
    Tres = Tuple{Int,String,String,String}
    seq = "GGGUUUAAAAAAACCCAAAUUUU"
    dbn = "(((...[[[[...)))...]]]]"
    for kwargs in [
        (; ),
        (; args=`-m`),
        ]
        res = run_RemovePseudoknots(seq, dbn; kwargs...)
        @test res isa Tres
    end
end

@testset "run_stochastic" begin
    Tres = Tuple{Int,String,String,String}
    seq = "GGGAAAACCC"

    for kwargs in [
        (; ),
        (; args=`-s 42`),
        ]
        res = run_stochastic(seq; kwargs...)
        @test res isa Tres
    end
end
