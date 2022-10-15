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
energy(seq::AbstractString, dbn::AbstractString; cmdline_opts=String[]) =
    first(energy(seq, [dbn]; cmdline_opts))
function energy(seq::AbstractString, dbns::Vector{<:AbstractString};
                cmdline_opts=String[])
    exitcode, res, out, err = efn2(seq, dbns; cmdline_opts)
    T = typeof(0.0 * UNIT_EN)
    energies = Tuple{T,T}[]
    try
        if exitcode != 0
            error("efn2 program returned nonzero exit status")
        end
        for line in eachline(IOBuffer(res))
            a = split(line)
            if length(a) != 7 || a[1] != "Structure:" || a[3] != "Energy" || a[4] != "=" || a[6] != "±"
                error("error parsing result line: $line")
            end
            en = en_stddev = 0.0
            try
                en = parse(Float64, a[5])
                en_stddev = parse(Float64, a[7])
            catch
                println("error parsing result line: $line")
                rethrow()
            end
            push!(energies, (en * UNIT_EN, en_stddev * UNIT_EN))
        end
        if length(energies) == 0
            error("no energies parsed")
        end
    catch
        println("contents of resultfile of efn2:")
        println(res, "\n")
        println("stdout of efn2:")
        println(out, "\n")
        println("stderr of efn2:")
        println(err, "\n")
        rethrow()
    end
    return energies
end

"""
    efn2(seq, dbn; cmdline_opts=String[]) -> exitcode, res, out, err
    efn2(seq, dbns; cmdline_opts=String[])

Run the `efn2` program from RNAstructure. Returns the exitcode of the
`efn2` program, the contents of the results file `res`, the stdout
`out`, and stderr `err` output as Strings.
"""
efn2(seq::AbstractString, dbn::AbstractString; cmdline_opts=String[]) =
    efn2(seq, [dbn]; cmdline_opts)
function efn2(seq::AbstractString,
              dbns::Vector{<:AbstractString};
              cmdline_opts=String[])
    exitcode = 0
    res = out = err = ""
    mktemp() do dbnpath, _
        mktemp() do respath, _
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
            buf_out = IOBuffer()
            buf_err = IOBuffer()
            cmd = `$(RNAstructure_jll.efn2()) $dbnpath $respath $cmdline_opts`
            r = run(pipeline(ignorestatus(cmd); stdout=buf_out, stderr=buf_err))
            exitcode = r.exitcode
            out = String(take!(buf_out))
            err = String(take!(buf_err))
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

end # module RNAstructure
