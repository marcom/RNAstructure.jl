module RNAstructure

import RNAstructure_jll
using Unitful: @u_str
import FoldRNA

export energy

include("ct-format.jl")

const UNIT_EN = u"kcal/mol"

function __init__()
    ENV["DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "data_tables")
end

"""
    energy(seq, dbn; cmdline_opts=String[])

Calculate free energy of folding of a nucleic acid sequence `seq`
folded into a secondary structure `dbn` given in dot-bracket notation
with the `efn2` program from RNAstructure.  Additional command-line
arguments can be passed with `cmdline_opts`.

Returns the free energy of folding and uncertainty.
"""
function energy(seq::AbstractString, dbn::AbstractString;
                cmdline_opts::Vector{<:Any}=String[])
    # TODO: better logic/error/exception handling and cleanup
    ctpath, _ = mktemp()
    respath, _ = mktemp()

    open(ctpath, "w") do io
        print_ct_format(io, dbn; seq)
    end

    # print(read(ctpath, String))

    buf_out = IOBuffer()
    buf_err = IOBuffer()
    run(pipeline(`$(RNAstructure_jll.efn2()) $ctpath $respath $cmdline_opts`;
                 stdout=buf_out, stderr=buf_err))
    out = String(take!(buf_out))
    err = String(take!(buf_err))
    res = read(respath, String)
    a = split(res)
    if length(a) != 7 || a[1] != "Structure:" || a[3] != "Energy" || a[4] != "=" || a[6] != "±"
        println("stdout of efn2:")
        println(out, "\n")
        println("stderr of efn2:")
        println(err, "\n")
        error("error parsing result: $res")
    end
    en = en_stddev = 0.0
    try
        en = parse(Float64, a[5])
        en_stddev = parse(Float64, a[7])
    catch e
        println("error parsing result: $res")
        throw(e)
    end

    rm(respath)
    rm(ctpath)
    return en * UNIT_EN, en_stddev * UNIT_EN
end

end # module RNAstructure
