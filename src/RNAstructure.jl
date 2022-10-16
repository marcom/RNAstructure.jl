module RNAstructure

import RNAstructure_jll
using Unitful: @u_str

export energy, design

const UNIT_EN = u"kcal/mol"

function __init__()
    ENV["DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "data_tables")
end

function _runcmd(cmd::Cmd)
    buf_out = IOBuffer()
    buf_err = IOBuffer()
    r = run(pipeline(ignorestatus(cmd); stdout=buf_out, stderr=buf_err))
    exitcode = r.exitcode
    out = String(take!(buf_out))
    err = String(take!(buf_err))
    return exitcode, out, err
end

_write_dbn_fasta(path::AbstractString, seq::AbstractString, dbn::AbstractString) =
    _write_dbn_fasta(path, seq, [dbn])

_write_dbn_fasta(path::AbstractString, seq::AbstractString) =
    _write_dbn_fasta(path, seq, String[])

function _write_dbn_fasta(path::AbstractString, seq::AbstractString,
                          dbns::Vector{<:AbstractString})
    # print dbn file format
    # dbn format:
    # >title
    # SEQUENCE
    # STRUCTURE1
    # STRUCTURE2...
    open(path, "w") do io
        println(io, ">")
        println(io, seq)
        for dbn in dbns
            println(io, dbn)
        end
    end
end

"""
    energy(seq, dbn; [verbose, cmdargs]) -> energy, uncertainty
    energy(seq, dbns; [verbose, cmdargs]) -> [energy, uncertainty]

Calculate free energy of folding of a nucleic acid sequence `seq`
folded into a secondary structure `dbn` (or an array of secondary
structures `dbns`) given in dot-bracket notation with the `efn2`
program from RNAstructure.  Additional command-line arguments can be
passed with `cmdargs`.

Returns the free energy of folding and experimental uncertainty.

See the [RNAstructure efn2
manpage](https://rna.urmc.rochester.edu/Text/efn2.html) for details on
command-line arguments that can be passed as `cmdargs`.
"""
function energy(seq::AbstractString, dbn::AbstractString;
                verbose::Bool=false, cmdargs=``)
    return first(energy(seq, [dbn]; verbose, cmdargs))
end

function energy(seq::AbstractString, dbns::Vector{<:AbstractString};
                verbose::Bool=false, cmdargs=``)
    exitcode, res, out, err = efn2(seq, dbns; cmdargs)
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
    if verbose
        println("stdout of efn2:")
        println(out, "\n")
        println("stderr of efn2:")
        println(err, "\n")
    end
    return energies
end

"""
    efn2(seq, dbn; [cmdargs]) -> exitcode, res, out, err
    efn2(seq, dbns; [cmdargs])

Run the `efn2` program from RNAstructure. Returns the exitcode of the
`efn2` program, the contents of the results file `res`, the stdout
`out`, and stderr `err` output as Strings.

See the [RNAstructure efn2
manpage](https://rna.urmc.rochester.edu/Text/efn2.html) for details on
command-line arguments that can be passed as `cmdargs`.
"""
efn2(seq::AbstractString, dbn::AbstractString; cmdargs=``) =
    efn2(seq, [dbn]; cmdargs)

function efn2(seq::AbstractString,
              dbns::Vector{<:AbstractString};
              cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do dbnpath, _
        mktemp() do respath, _
            _write_dbn_fasta(dbnpath, seq, dbns)
            cmd = `$(RNAstructure_jll.efn2()) $dbnpath $respath $cmdargs`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    design(target_dbn; [verbose, cmdargs]) -> seq, seed

Design sequences that will fold into the secondary structure
`target_dbn` given in dot-bracket notation.

See the [RNAstructure design
manpage](https://rna.urmc.rochester.edu/Text/design.html) for details
on command-line arguments that can be passed as `cmdargs`.
"""
function design(target_dbn::AbstractString;
                verbose::Bool=false, cmdargs=``)
    exitcode = 0
    out = err = ""
    mktemp() do dbnpath, _
        _write_dbn_fasta(dbnpath, "N"^length(target_dbn), target_dbn)
        cmd = `$(RNAstructure_jll.design()) $dbnpath $cmdargs`
        exitcode, out, err = _runcmd(cmd)
    end
    if verbose || exitcode != 0
        println("stdout of design:")
        println(out, "\n")
        println("stderr of design:")
        println(err, "\n")
    end
    if exitcode != 0
        error("design returned non-zero exit status")
    end
    seq = match(r"\nResult= (\S+)\n", out).captures[1] |> String
    seed = match(r"\sRandomSeed:\s+(\S+)\s", out).captures[1] |> String
    return (; seq, seed)
end

"""
    fold(seq; [verbose, cmdargs]) -> res, out, err

Run the `Fold` program from RNAstructure.

See the [RNAstructure Fold
manpage](https://rna.urmc.rochester.edu/Text/Fold.html) for details on
command-line arguments that can be passed as `cmdargs`.
"""
function fold(seq::AbstractString; verbose::Bool=false, cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do dbnpath, _
            _write_dbn_fasta(dbnpath, seq)
            cmd = `$(RNAstructure_jll.Fold()) $dbnpath $respath $cmdargs`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    if verbose || exitcode != 0
        println("stdout of Fold:")
        println(out, "\n")
        println("stderr of Fold:")
        println(err, "\n")
    end
    if exitcode != 0
        error("Fold returned non-zero exit status")
    end
    return res, out, err
end

end # module RNAstructure
