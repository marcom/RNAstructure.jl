module RNAstructure

import RNAstructure_jll
using Unitful: @u_str

export energy

const UNIT_EN = u"kcal/mol"

function __init__()
    ENV["DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "data_tables")
end

"""
    energy(seq, dbn; cmdline_opts=String[])
    energy(seq, dbns; cmdline_opts=String[])

Calculate free energy of folding of a nucleic acid sequence `seq`
folded into a secondary structure `dbn` (or an array of secondary
structures `dbns`) given in dot-bracket notation with the `efn2`
program from RNAstructure.  Additional command-line arguments can be
passed with `cmdline_opts`.

Returns the free energy of folding and experimental uncertainty.
"""
energy(seq::AbstractString, dbn::AbstractString; kwargs...) =
    first(energy(seq, [dbn]; kwargs...))
function energy(seq::AbstractString,
                dbns::Vector{<:AbstractString};
                cmdline_opts::Vector{<:Any}=String[])
    # TODO: better logic/error/exception handling and cleanup
    dbnpath, _ = mktemp()
    respath, _ = mktemp()

    # print dbn file format
    # dbn format:
    # >title
    # SEQUENCE
    # STRUCTURE1
    # STRUCTURE2...
    open(dbnpath, "w") do io
        println(io, ">")
        println(io, seq)
        for dbn in dbns
            println(io, dbn)
        end
    end

    out = err = ""
    try
        buf_out = IOBuffer()
        buf_err = IOBuffer()
        run(pipeline(`$(RNAstructure_jll.efn2()) $dbnpath $respath $cmdline_opts`;
                     stdout=buf_out, stderr=buf_err))
        out = String(take!(buf_out))
        err = String(take!(buf_err))

        T = typeof(0.0 * UNIT_EN)
        energies = Tuple{T,T}[]
        for line in eachline(respath)
            a = split(line)
            if length(a) != 7 || a[1] != "Structure:" || a[3] != "Energy" || a[4] != "=" || a[6] != "±"
                error("error parsing result line: $line")
            end
            en = en_stddev = 0.0
            try
                en = parse(Float64, a[5])
                en_stddev = parse(Float64, a[7])
            catch e
                println("error parsing result line: $line")
                rethrow(e)
            end
            push!(energies, (en * UNIT_EN, en_stddev * UNIT_EN))
        end
        if length(energies) == 0
            error("no energies parsed")
        end
    catch e
        println("stdout of efn2:")
        println(out, "\n")
        println("stderr of efn2:")
        println(err, "\n")
        rethrow(e)
    end

    rm(respath)
    rm(dbnpath)
    return energies
end

end # module RNAstructure
