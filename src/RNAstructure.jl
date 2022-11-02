module RNAstructure

import RNAstructure_jll
using Unitful: Quantity, @u_str, uconvert, ustrip

export bpp, dbn2ct, design, energy, ensemble_defect, mea, mfe, partfn,
    prob_of_structure, remove_pknots, sample_structures, subopt

const UNIT_EN = u"kcal/mol"

include("ct-format.jl")
include("pairtable-to-dbn.jl")

function __init__()
    if !haskey(ENV, "DATAPATH")
        ENV["DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "data_tables")
    else
        @info "RNAstructure: energy params already set, DATAPATH=$(ENV["DATAPATH"])"
    end
    # TODO: set OMP_NUM_THREADS env var for smp programs (number of threads to use)
    # TODO: set CYCLEFOLD_DATAPATH env var to
    #       RNAstructure/CycleFold/datafiles path (is this in R2R_jll?)
    return nothing
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

function _parse_bpp_file!(p, bpp_str::AbstractString)
    # basepair prob text file format (from `ProbabilityPlot --text`)
    # <n>
    # i           j	       -log10(Probability)
    # <idx_i>	  <idx_j>      <-log10_prob>
    # ...
    re_emptyline = r"^\s*$"
    io = IOBuffer(bpp_str)
    want_firstline = true
    want_secondline = false
    for line in eachline(io)
        if occursin(re_emptyline, line)
            continue
        end
        if want_firstline
            n = parse(Int, line)
            want_firstline = false
            want_secondline = true
            continue
        end
        if want_secondline
            if split(line) != ["i", "j", "-log10(Probability)"]
                error("expected header line, got: $line")
            end
            want_secondline = false
            continue
        end
        i_str, j_str, m_log10_prob_str = split(line)
        i = parse(Int, i_str)
        j = parse(Int, j_str)
        pij = 10.0^(- parse(Float64, m_log10_prob_str))
        p[i,j] = pij
        p[j,i] = pij
    end
end

"""
    bpp(seq; [verbose, cmdargs]) -> basepair_prob

Calculate basepair probabilities for an RNA sequence `seq`.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function bpp(seq::AbstractString;
             verbose::Bool=false, cmdargs=``)
    # TODO: min and max value to save
    # TODO: return type as sparse matrix? (see LinearFold)
    n = length(seq)
    pij = zeros(n, n)
    mktemp() do pf_savefile, _
        exitcode, out, err = run_partition!(pf_savefile, seq; cmdargs)
        want_help = ("-h" in cmdargs) || ("--help" in cmdargs)
        if verbose || exitcode != 0 || want_help
            println("stdout of partition:")
            println(out)
            println("stderr of partition:")
            println(err)
        end
        if exitcode != 0
            error("partition returned non-zero exit status ($exitcode)")
        end
        if want_help
            # TODO: we throw a error here, because most RNAstructure
            # programs return a non-zero exit status when help is
            # requested, but `partition` doesn't
            error("help string requested")
        end
        exitcode, res, out, err = run_ProbabilityPlot(pf_savefile; cmdargs=`--text`)
        if verbose || exitcode != 0
            println("stdout of ProbabilityPlot:")
            println(out)
            println("stderr of ProbabilityPlot:")
            println(err)
        end
        if exitcode != 0
            error("ProbabilityPlot returned non-zero exit status ($exitcode)")
        end
        _parse_bpp_file!(pij, res)
    end
    return pij
end

"""
    dbn2ct(dbn; [verbose]) -> ct::String
    dbn2ct(seq, dbn; [verbose]) -> ct::String

Convert secondary structure in dot-bracket format `dbn` to ct format.
If no sequence is given, the sequence will be all 'N'.
"""
function dbn2ct(seq::AbstractString, dbn::AbstractString; verbose::Bool=false)
    exitcode, ct, out, err = run_dot2ct(seq, dbn)
    if verbose || exitcode != 0
        println("stdout of dot2ct:")
        println(out)
        println("stderr of dot2ct:")
        println(err)
    end
    exitcode == 0 || error("dot2ct returned non-zero exit status ($exitcode)")
    return ct
end

dbn2ct(dbn::AbstractString; verbose::Bool=false) =
    dbn2ct("N"^length(dbn), dbn; verbose)

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
            # lines are of the form
            # Structure: 1   Energy = -1.2 ± 0.1
            # or
            # Structure: 1   Energy = -0.2
            a = split(line)
            if length(a) == 5
                # if a[1] != "Structure:" || a[3] != "Energy" || a[4] != "="
                #     error("wrong form for result line: $line")
                # end
                error("result lines of this form not yet supported: $line")
            elseif length(a) == 7
                if a[1] != "Structure:" || a[3] != "Energy" || a[4] != "=" || a[6] != "±"
                    error("wrong form for result line: $line")
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
            else
                error("wrong number of tokens on result line: $line")
            end
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
    mea(seq; [verbose, cmdargs_partition, cmdargs_maxexpect]) -> mea_structures

Calculate maximum expected accuracy (MEA) structures for an RNA
sequence `seq`. Returns secondary structures in dot-bracket notation
as a vector of Strings.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html) for
details on command-line arguments that can be passed as
`cmdargs_partition`.

See the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html) for
details on command-line arguments that can be passed as
`cmdargs_maxexpect`.
"""
function mea(seq::AbstractString;
             verbose::Bool=false, cmdargs_partition=``, cmdargs_maxexpect=``)
    want_help = ("-h" in cmdargs_partition) || ("-h" in cmdargs_maxexpect) ||
        ("--help" in cmdargs_partition) || ("--help" in cmdargs_maxexpect)
    if want_help
        exitcode, out, err = run_partition!("", ""; cmdargs=`-h`)
        println("partition\n", out, err)
        exitcode, out, err = run_MaxExpect(""; cmdargs=`-h`)
        println("MaxExpect\n", out, err)
        # TODO: we throw a error here, because most RNAstructure
        # programs return a non-zero exit status when help is
        # requested, but `partition` doesn't
        error("help string requested")
    end
    res = ""
    mktemp() do pf_savefile, _
        exitcode, out, err = run_partition!(pf_savefile, seq; cmdargs=cmdargs_partition)
        if verbose || exitcode != 0
            println("stdout of partition:")
            println(out)
            println("stderr of partition:")
            println(err)
        end
        exitcode == 0 || error("partition returned non-zero exit status ($exitcode)")
        exitcode, res, out, err = run_MaxExpect(pf_savefile; cmdargs=cmdargs_maxexpect)
        if verbose || exitcode != 0
            println("stdout of MaxExpect:")
            println(out)
            println("stderr of MaxExpect:")
            println(err)
        end
        exitcode == 0 || error("MaxExpect returned non-zero exit status ($exitcode)")
    end
    ct_structs = parse_ct_format(res)
    dbns = [pairtable_to_dbn(pt) for (_, _, pt) in ct_structs]
    return dbns
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
    partfn(seq; [verbose, cmdargs]) -> ensemble_energy

Calculate the partition function for RNA sequence `seq` and return the
ensemble energy.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function partfn(seq; verbose::Bool=false, cmdargs=``)
    # TODO: if partition function save files are used, i think a lot
    #       more `cmdargs` would be possible (would mean using both
    #       run_partition and then run_EnsembleEnergy)
    exitcode, out, err = run_EnsembleEnergy(seq; cmdargs)
    m = match(r"\nEnsemble energy.*:(.*) kcal/mol\n", out)
    if verbose || exitcode != 0 || isnothing(m)
        println("stdout of EnsembleEnergy:")
        println(out, "\n")
        println("stderr of EnsembleEnergy:")
        println(err, "\n")
    end
    if exitcode != 0
        error("EnsembleEnergy returned non-zero exit status")
    end
    if isnothing(m)
        error("couldn't find ensemble energy in output")
    end
    en = parse(Float64, m.captures[1])
    return en * u"kcal/mol"
end

"""
    prob_of_structure(seq, dbn; [cmdargs])

Calculates the probability of a given RNA sequence `seq` folding into
a secondary structure `dbn` given in dot-bracket notation.

Note: currently the temperature is fixed at 37°C and cannot be
changed.

The supported `cmdargs` are those common to `energy` and `partfn`.
"""
function prob_of_structure(seq::AbstractString, dbn::AbstractString;
                           cmdargs=``)
    # TODO
    # - make temperature a kwarg once partfn supports it
    # - support cmdargs once partfn uses partition (common kwargs of
    #   partition and efn2)
    temperature = uconvert(u"K", 37u"°C")
    # t = float(ustrip(temperature))
    # cmdargs = `$cmdargs -T $t`
    en, _ = energy(seq, dbn; cmdargs)
    pf = partfn(seq; cmdargs)  # pf == -RT * log(Q)
    R = uconvert(u"kcal/K/mol", 1.98720425864083u"cal/K/mol")
    RT = R * temperature
    log_p = (-en + pf) / RT
    return exp(log_p)
end

"""
    remove_pknots(dbn; [verbose]) -> dbn

Remove pseudoknots from a secondary structure `dbn` in dot-bracket
format by removing the fewest possible basepairs.
"""
function remove_pknots(dbn::AbstractString; verbose::Bool=false)
    # TODO: multiple structures at a time possible?
    # TODO: return pknot-free struct with lowest mfe (don't use -m,
    #       accept more cmdargs?)
    exitcode, res, out, err = run_RemovePseudoknots("N"^length(dbn), dbn;
                                                    cmdargs=`-m`)
    if verbose || exitcode != 0
        println("stdout of RemovePseudoknots:")
        println(out)
        println("stderr of RemovePseudoknots:")
        println(err)
    end
    exitcode == 0 || error("RemovePseudoknots returned non-zero exit status ($exitcode)")
    ct_structs = parse_ct_format(res)
    dbns = [pairtable_to_dbn(pt) for (_, _, pt) in ct_structs]
    return first(dbns)
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
    subopt(seq; [verbose, cmdargs]) -> subopt_en_structures

Calculate suboptimal structures for an RNA sequence `seq`. Returns a
vector of secondary structures in dot-bracket notation and their
energies.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for
details on command-line arguments that can be passed as `cmdargs`.

In particular, the `--maximum`, `--percent`, and `--window` arguments
can be used to increase or decrease the number of secondary structures
returned.
"""
function subopt(seq::AbstractString; verbose::Bool=false, cmdargs=``)
    if ("-h" in cmdargs) || ("--help" in cmdargs)
        exitcode, res, out, err = run_Fold(""; cmdargs=`-h`)
        println(out, err)
        error("help string requested")
    end
    exitcode, res, out, err = run_Fold(seq; cmdargs)
    if verbose || exitcode != 0
        println("stdout of Fold:")
        println(out)
        println("stderr of Fold:")
        println(err)
    end
    # TODO: make energy parsing a global helper fn
    get_energy(str) = try parse(Float64, split(str, "=")[2]) catch; 0.0 end * u"kcal/mol"
    exitcode == 0 || error("Fold returned non-zero exit status ($exitcode)")
    ct_structs = parse_ct_format(res)
    en_dbns = [(pairtable_to_dbn(pt), get_energy(title)) for (title, _, pt) in ct_structs]
    return en_dbns
end

"""
    run_dot2ct(seq, dbn; [cmdargs]) -> exitcode, res, out, err
    run_dot2ct(dbn; [cmdargs]) -> exitcode, res, out, err

Run the `dot2ct` program from RNAstructure.

See the [RNAstructure dot2ct
documentation](https://rna.urmc.rochester.edu/Text/dot2ct.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function run_dot2ct(seq::AbstractString, dbn::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do dbnpath, _
            _write_dbn_fasta(dbnpath, seq, dbn)
            cmd = `$(RNAstructure_jll.dot2ct()) $dbnpath $respath $cmdargs`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

run_dot2ct(dbn::AbstractString; cmdargs=``) =
    run_dot2ct("N"^length(dbn), dbn; cmdargs)

"""
    run_draw(dbn, [seq]; [cmdargs]) -> exitcode, res, out, err

Run the `draw` program from RNAstructure.

See the [RNAstructure draw
documentation](https://rna.urmc.rochester.edu/Text/draw.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function run_draw(dbn::AbstractString, seq::AbstractString;
                  cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do dbnpath, _
            _write_dbn_fasta(dbnpath, seq, dbn)
            cmd = `$(RNAstructure_jll.draw()) $dbnpath $respath $cmdargs`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
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
    run_EnsembleEnergy(seq; [cmdargs]) -> exitcode, out, err

Run the `EnsembleEnergy` program from RNAstructure. Returns the
exitcode of the `EnsembleEnergy` program, the stdout `out`, and stderr
`err` output as Strings.

See the [RNAstructure EnsembleEnergy
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function run_EnsembleEnergy(seq::AbstractString; cmdargs=``)
    exitcode = 0
    out = err = ""
    mktemp() do seqpath, _
        _write_dbn_fasta(seqpath, seq)
        cmd = `$(RNAstructure_jll.EnsembleEnergy()) $seqpath --sequence $cmdargs`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
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
    run_MaxExpect(seq; [cmdargs]) -> exitcode, res, out, err

Run the `MaxExpect` program from RNAstructure.

See the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html) for
details on command-line arguments that can be passed as `cmdargs`.
"""
function run_MaxExpect(pf_savefile::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        cmd = `$(RNAstructure_jll.MaxExpect()) $pf_savefile $respath $cmdargs`
        exitcode, out, err = _runcmd(cmd)
        res = read(respath, String)
    end
    return exitcode, res, out, err
end

"""
    run_partition!(pf_savefile, seq; [cmdargs]) -> exitcode, out, err

Run the `partition` program from RNAstructure. The partition function
save file is saved to the path given by `pf_savefile`.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function run_partition!(pf_savefile::AbstractString, seq::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do seqpath, _
        _write_dbn_fasta(seqpath, seq)
        cmd = `$(RNAstructure_jll.partition()) $seqpath $pf_savefile $cmdargs`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
end

"""
    run_ProbabilityPlot(pf_savefile; [cmdargs]) -> exitcode, res, out, err

Run the `ProbabilityPlot` program from RNAstructure.  The partition
function save file `pf_savefile` is the path to a file previously
generated with `run_partition!`.

See the [RNAstructure ProbabilityPlot
documentation](https://rna.urmc.rochester.edu/Text/ProbabilityPlot.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function run_ProbabilityPlot(pf_savefile::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        cmd = `$(RNAstructure_jll.ProbabilityPlot()) $pf_savefile $respath $cmdargs`
        exitcode, out, err = _runcmd(cmd)
        res = read(respath, String)
    end
    return exitcode, res, out, err
end

"""
    run_RemovePseudoknots(seq, dbn; [cmdargs]) -> exitcode, res, out, err

Run the `RemovePseudoknots` program from RNAstructure.

See the [RNAstructure RemovePseudoknots
documentation](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html)
for details on command-line arguments that can be passed as `cmdargs`.
"""
function run_RemovePseudoknots(seq::AbstractString, dbn::AbstractString; cmdargs=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do ctpath, _
            ps, ct, _... = run_dot2ct(seq, dbn)
            ps == 0 || error("dot2ct returned non-zero exit status")
            open(ctpath, "w") do io
                write(io, ct)
            end
            cmd = `$(RNAstructure_jll.RemovePseudoknots()) $ctpath $respath $cmdargs`
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
