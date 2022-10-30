module RNAstructure

import RNAstructure_jll
using Unitful: @u_str

export design, energy, ensemble_defect, mfe, sample_structures

const UNIT_EN = u"kcal/mol"

include("parse-ct-format.jl")
include("pairtable-to-dbn.jl")

function __init__()
    ENV["DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "data_tables")
    # TODO: set OMP_NUM_THREADS env var for smp programs (number of threads to use)
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
    design(target_dbn; [verbose, cmdargs]) -> seq, seed

Design sequences that will fold into the secondary structure
`target_dbn` given in dot-bracket notation.

See the [RNAstructure design
documentation](https://rna.urmc.rochester.edu/Text/design.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function design(target_dbn::AbstractString;
                verbose::Bool=false, cmdargs=``)
    # TODO: split this function into `run_design` and `design`
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
    energy(seq, dbn; [verbose, cmdargs]) -> energy, uncertainty
    energy(seq, dbns; [verbose, cmdargs]) -> [(energy, uncertainty), ...]

Calculate free energy of folding of a nucleic acid sequence `seq`
folded into a secondary structure `dbn` (or an array of secondary
structures `dbns`) given in dot-bracket notation with the `efn2`
program from RNAstructure.  Additional command-line arguments can be
passed with `cmdargs`.

Returns the free energy of folding and experimental uncertainty.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function energy(seq::AbstractString, dbn::AbstractString;
                verbose::Bool=false, cmdargs=``)
    return first(energy(seq, [dbn]; verbose, cmdargs))
end

function energy(seq::AbstractString, dbns::Vector{<:AbstractString};
                verbose::Bool=false, cmdargs=``)
    exitcode, res, out, err = run_efn2(seq, dbns; cmdargs)
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
    ensemble_defect(seq, dbn; [verbose, cmdargs]) -> ed, ned
    ensemble_defect(seq, dbns; [verbose, cmdargs]) -> [(ed, ned), ...]

Calculates the ensemble defect for a sequence `seq` and a secondary
structure `dbn` or multiple secondary structures `dbns` given in
dot-bracket notation.  Additional command-line arguments can be passed
with `cmdargs` to the `EDcalculator` program from RNAstructure.

Returns the ensemble defect and normalised ensemble defect.

See the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function ensemble_defect(seq::AbstractString, dbn::AbstractString;
                         verbose::Bool=false, cmdargs=``)
    return first(ensemble_defect(seq, [dbn]; verbose, cmdargs))
end

function ensemble_defect(seq::AbstractString, dbns::Vector{<:AbstractString};
                         verbose::Bool=false, cmdargs=``)
    exitcode, out, err = run_EDcalculator(seq, dbns; cmdargs)
    eds = Tuple{Float64,Float64}[]
    # Output lines have this form:
    # Structure 1: Ensemble_Defect =	3.17124		Normalized_ED =	0.352361
    try
        regex = r"^Structure\s+\d+:\s*Ensemble_Defect\s*=\s*(\S+)\s+Normalized_ED\s*=\s*(\S+)\s*$"m
        for m in eachmatch(regex, out)
            str_ed, str_ned = m.captures
            ed = parse(Float64, str_ed)
            ned = parse(Float64, str_ned)
            push!(eds, (ed, ned))
        end
    catch
        println("stdout of EDcalculator:")
        println(out, "\n")
        println("stderr of EDcalculator:")
        println(err, "\n")
        rethrow()
    end
    if verbose || exitcode != 0
        println("stdout of EDcalculator:")
        println(out, "\n")
        println("stderr of EDcalculator:")
        println(err, "\n")
    end
    if exitcode != 0
        error("EDcalculator returned non-zero exit status")
    end
    return eds
end

"""
    mfe(seq; [verbose, cmdargs]) -> energy, mfe_structure

Calculate the minimum free energy (MFE) structure of an RNA sequence
`seq` by calling the `Fold` program from RNAstructure.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function mfe(seq; verbose::Bool=false, cmdargs=``)
    exitcode, res, out, err = run_Fold(seq; cmdargs=`-mfe $cmdargs`)
    if verbose || exitcode != 0
        println("result file of Fold:")
        println(res, "\n")
        println("stdout of Fold:")
        println(out, "\n")
        println("stderr of Fold:")
        println(err, "\n")
    end
    if exitcode != 0
        error("Fold returned non-zero exit status")
    end
    ct_structs = parse_ct_format(res)
    if length(ct_structs) != 1
        error("expected exactly one structure, got $(length(ct_structs))\n",
              "structures are: $ct_structs")
    end
    title, seq, pairtable = ct_structs[1]
    a = split(title, "=")
    en = if length(a) == 2
        parse(Float64, a[2]) * u"kcal/mol"
    elseif length(a) == 1
        0.0u"kcal/mol"
    else
        error("could not parse title line to find energy: $title")
    end
    return en, pairtable_to_dbn(pairtable)
end

"""
    sample_structures(seq; [verbose, cmdargs]) -> dbn_structures::Vector{String}

Sample secondary structures from the Boltzmann ensemble of structures
for the RNA sequence `seq`. Returns an array of secondary structures
in dot-bracket notation.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function sample_structures(seq; verbose::Bool=false, cmdargs=``)
    exitcode, res, out, err = run_stochastic(seq; cmdargs=`$cmdargs`)
    if verbose || exitcode != 0
        println("result file of stochastic:")
        println(res, "\n")
        println("stdout of stochastic:")
        println(out, "\n")
        println("stderr of stochastic:")
        println(err, "\n")
    end
    if exitcode != 0
        error("stochastic returned non-zero exit status")
    end
    ct_structs = parse_ct_format(res)
    dbns = [pairtable_to_dbn(pt) for (_, _, pt) in ct_structs]
    return dbns
end

"""
    run_EDcalculator(seq, dbn; [cmdargs]) -> exitcode, out, err
    run_EDcalculator(seq, dbns; [cmdargs])

Run the `EDcalculator` program from RNAstructure.

See the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function run_EDcalculator(seq::AbstractString, dbn::AbstractString;
                          cmdargs=``)
    return run_EDcalculator(seq, [dbn]; cmdargs)
end

function run_EDcalculator(seq::AbstractString, dbns::Vector{<:AbstractString};
                          cmdargs=``)
    exitcode = 0
    out = err = ""
    mktemp() do dbnpath, _
        _write_dbn_fasta(dbnpath, seq, dbns)
        cmd = `$(RNAstructure_jll.EDcalculator()) $dbnpath $cmdargs`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
end

"""
    run_efn2(seq, dbn; [cmdargs]) -> exitcode, res, out, err
    run_efn2(seq, dbns; [cmdargs])

Run the `efn2` program from RNAstructure. Returns the exitcode of the
`efn2` program, the contents of the results file `res`, the stdout
`out`, and stderr `err` output as Strings.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
run_efn2(seq::AbstractString, dbn::AbstractString; cmdargs=``) =
    run_efn2(seq, [dbn]; cmdargs)

function run_efn2(seq::AbstractString,
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
    run_Fold(seq; [cmdargs]) -> exitcode, res, out, err

Run the `Fold` program from RNAstructure.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function run_Fold(seq::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do seqpath, _
            _write_dbn_fasta(seqpath, seq)
            cmd = `$(RNAstructure_jll.Fold()) $seqpath $respath $cmdargs`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_stochastic(seq; [cmdargs]) -> exitcode, res, out, err

Run the `stochastic` program from RNAstructure.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function run_stochastic(seq::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do seqpath, _
            _write_dbn_fasta(seqpath, seq)
            cmd = `$(RNAstructure_jll.stochastic()) $seqpath $respath --sequence $cmdargs`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

end # module RNAstructure
