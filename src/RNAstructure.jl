module RNAstructure

import RNAstructure_jll
using Unitful: Quantity, @u_str, uconvert, ustrip

export bpp, ct2dbn, dbn2ct, design, energy, ensemble_defect, mea, mfe,
    partfn, plot, prob_of_structure, remove_pknots, sample_structures,
    subopt, subopt_all

const UNIT_EN = u"kcal/mol"

include("ct-format.jl")
include("pairtable-to-dbn.jl")
include("plot.jl")

function __init__()
    if !haskey(ENV, "DATAPATH")
        ENV["DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "data_tables")
    else
        @info "RNAstructure: energy params already set, DATAPATH=$(ENV["DATAPATH"])"
    end
    if !haskey(ENV, "CYCLEFOLD_DATAPATH")
        ENV["CYCLEFOLD_DATAPATH"] = joinpath(RNAstructure_jll.artifact_dir, "CycleFold", "datafiles")
    else
        @info "RNAstructure: CycleFold energy params already set, CYCLEFOLD_DATAPATH=$(ENV["CYCLEFOLD_DATAPATH"])"
    end
    # TODO: set OMP_NUM_THREADS env var for smp programs (number of threads to use)
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
    # print "fasta-dbn" file format
    # format:
    # >title        (on one line)
    # SEQUENCE      (on one line)
    # STRUCTURE1    (on one line)
    # STRUCTURE<N>...
    open(path, "w") do io
        println(io, ">")
        println(io, seq)
        for dbn in dbns
            println(io, dbn)
        end
    end
end

# returns: title, seq, [dbns]
function _parse_dbn_fasta(dbnfasta::AbstractString)
    lines = readlines(IOBuffer(dbnfasta))
    if length(lines) < 3
        throw(ArgumentError("expected at least 3 input lines, only found $(length(lines))"))
    end
    title_arr = collect(lines[1])
    seq = lines[2]
    dbns = lines[3:end]
    if title_arr[1] != '>'
        throw(ArgumentError("expected title line to start with '>'"))
    end
    title = join(title_arr[2:end])
    return title, seq, dbns
end

function _write_fasta(path::AbstractString, title_seq_iterable)
    # print (multi) fasta file format
    # >title1
    # SEQUENCE1
    # >title2
    # SEQUENCE2
    # ...
    open(path, "w") do io
        for (title, seq) in title_seq_iterable
            println(io, ">", title, "\n", seq)
        end
    end
end

_write_fasta(path::AbstractString, seqs::Vector{<:AbstractString}) =
    _write_fasta(path, [("", s) for s in seqs])

_write_fasta(path::AbstractString, seq::AbstractString) =
    _write_fasta(path, [("", seq)])

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
    bpp(seq; [verbose, args]) -> basepair_prob

Calculate basepair probabilities for an RNA sequence `seq`.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html) for
details on command-line arguments that can be passed as `args`.
"""
function bpp(seq::AbstractString;
             verbose::Bool=false, args=``)
    # TODO: min and max value to save
    # TODO: return type as sparse matrix? (see LinearFold)
    n = length(seq)
    pij = zeros(n, n)
    mktemp() do pf_savefile, _
        exitcode, out, err = run_partition!(pf_savefile, seq; args)
        want_help = ("-h" in args) || ("--help" in args)
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
        exitcode, res, out, err = run_ProbabilityPlot(pf_savefile; args=`--text`)
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
    ct2dbn(ct::AbstractString, i::Integer=1; [verbose]) -> title, seq, [dbns]

Convert RNA secondary structure number `i` from `ct` format to
dot-bracket notation.
"""
function ct2dbn(ct::AbstractString, i::Integer=1; verbose::Bool=false)
    exitcode, res, out, err = run_ct2dot(ct, i)
    if verbose || exitcode != 0
        println("stdout of ct2dot:")
        println(out)
        println("stderr of ct2dot:")
        println(err)
    end
    exitcode == 0 || error("ct2dot returned non-zero exit status ($exitcode)")
    title, seq, dbns = _parse_dbn_fasta(res)
    return title, seq, dbns
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
    design(target_dbn; [verbose, args]) -> seq, seed

Design sequences that will fold into the secondary structure
`target_dbn` given in dot-bracket notation.

See the [RNAstructure design
documentation](https://rna.urmc.rochester.edu/Text/design.html) for
details on command-line arguments that can be passed as `args`.
"""
function design(target_dbn::AbstractString;
                verbose::Bool=false, args=``)
    # TODO: split this function into `run_design` and `design`
    exitcode = 0
    out = err = ""
    mktemp() do dbnpath, _
        _write_dbn_fasta(dbnpath, "N"^length(target_dbn), target_dbn)
        cmd = `$(RNAstructure_jll.design()) $dbnpath $args`
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
    energy(seq, dbn; [verbose, args]) -> energy, uncertainty
    energy(seq, dbns; [verbose, args]) -> [(energy, uncertainty), ...]

Calculate free energy of folding of a nucleic acid sequence `seq`
folded into a secondary structure `dbn` (or an array of secondary
structures `dbns`) given in dot-bracket notation with the `efn2`
program from RNAstructure.  Additional command-line arguments can be
passed with `args`.

Returns the free energy of folding and experimental uncertainty.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for
details on command-line arguments that can be passed as `args`.
"""
function energy(seq::AbstractString, dbn::AbstractString;
                verbose::Bool=false, args=``)
    return first(energy(seq, [dbn]; verbose, args))
end

function energy(seq::AbstractString, dbns::Vector{<:AbstractString};
                verbose::Bool=false, args=``)
    exitcode, res, out, err = run_efn2(seq, dbns; args)
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
    ensemble_defect(seq, dbn; [verbose, args]) -> ed, ned
    ensemble_defect(seq, dbns; [verbose, args]) -> [(ed, ned), ...]

Calculates the ensemble defect for a sequence `seq` and a secondary
structure `dbn` or multiple secondary structures `dbns` given in
dot-bracket notation.  Additional command-line arguments can be passed
with `args` to the `EDcalculator` program from RNAstructure.

Returns the ensemble defect and normalised ensemble defect.

See the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html)
for details on command-line arguments that can be passed as `args`.
"""
function ensemble_defect(seq::AbstractString, dbn::AbstractString;
                         verbose::Bool=false, args=``)
    return first(ensemble_defect(seq, [dbn]; verbose, args))
end

function ensemble_defect(seq::AbstractString, dbns::Vector{<:AbstractString};
                         verbose::Bool=false, args=``)
    exitcode, out, err = run_EDcalculator(seq, dbns; args)
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
    mea(seq; [verbose, args_partition, args_maxexpect]) -> mea_structures

Calculate maximum expected accuracy (MEA) structures for an RNA
sequence `seq`. Returns secondary structures in dot-bracket notation
as a vector of Strings.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html) for
details on command-line arguments that can be passed as
`args_partition`.

See the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html) for
details on command-line arguments that can be passed as
`args_maxexpect`.
"""
function mea(seq::AbstractString;
             verbose::Bool=false, args_partition=``, args_maxexpect=``)
    want_help = ("-h" in args_partition) || ("-h" in args_maxexpect) ||
        ("--help" in args_partition) || ("--help" in args_maxexpect)
    if want_help
        exitcode, out, err = run_partition!("", ""; args=`-h`)
        println("partition\n", out, err)
        exitcode, out, err = run_MaxExpect(""; args=`-h`)
        println("MaxExpect\n", out, err)
        # TODO: we throw a error here, because most RNAstructure
        # programs return a non-zero exit status when help is
        # requested, but `partition` doesn't
        error("help string requested")
    end
    res = ""
    mktemp() do pf_savefile, _
        exitcode, out, err = run_partition!(pf_savefile, seq; args=args_partition)
        if verbose || exitcode != 0
            println("stdout of partition:")
            println(out)
            println("stderr of partition:")
            println(err)
        end
        exitcode == 0 || error("partition returned non-zero exit status ($exitcode)")
        exitcode, res, out, err = run_MaxExpect(pf_savefile; args=args_maxexpect)
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
    mfe(seq; [verbose, args]) -> energy, mfe_structure

Calculate the minimum free energy (MFE) structure of an RNA sequence
`seq` by calling the `Fold` program from RNAstructure.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for
details on command-line arguments that can be passed as `args`.
"""
function mfe(seq; verbose::Bool=false, args=``)
    exitcode, res, out, err = run_Fold(seq; args=`-mfe $args`)
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
    partfn(seq; [verbose, args]) -> ensemble_energy

Calculate the partition function for RNA sequence `seq` and return the
ensemble energy.

See the [RNAstructure EnsembleEnergy
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html)
for details on command-line arguments that can be passed as `args`.
"""
function partfn(seq; verbose::Bool=false, args=``)
    # TODO: if partition function save files are used, i think a lot
    #       more `args` would be possible (would mean using both
    #       run_partition and then run_EnsembleEnergy)
    exitcode, out, err = run_EnsembleEnergy(seq; args)
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
    prob_of_structure(seq, dbn; [args])

Calculates the probability of a given RNA sequence `seq` folding into
a secondary structure `dbn` given in dot-bracket notation.

Note: currently the temperature is fixed at 37°C and cannot be
changed.

The supported `args` are those common to `energy` and `partfn`.
"""
function prob_of_structure(seq::AbstractString, dbn::AbstractString;
                           args=``)
    # TODO
    # - make temperature a kwarg once partfn supports it
    # - support args once partfn uses partition (common kwargs of
    #   partition and efn2)
    temperature = uconvert(u"K", 37u"°C")
    # t = float(ustrip(temperature))
    # args = `$args -T $t`
    en, _ = energy(seq, dbn; args)
    pf = partfn(seq; args)  # pf == -RT * log(Q)
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
    #       accept more args?)
    exitcode, res, out, err = run_RemovePseudoknots("N"^length(dbn), dbn;
                                                    args=`-m`)
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
    sample_structures(seq; [verbose, args]) -> dbn_structures::Vector{String}

Sample secondary structures from the Boltzmann ensemble of structures
for the RNA sequence `seq`. Returns an array of secondary structures
in dot-bracket notation.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html)
for details on command-line arguments that can be passed as `args`.
"""
function sample_structures(seq; verbose::Bool=false, args=``)
    exitcode, res, out, err = run_stochastic(seq; args=`$args`)
    want_help = ("-h" in args) || ("--help" in args)
    if verbose || exitcode != 0 || want_help
        println("result file of stochastic:")
        println(res, "\n")
        println("stdout of stochastic:")
        println(out, "\n")
        println("stderr of stochastic:")
        println(err, "\n")
    end
    exitcode == 0 || error("stochastic returned non-zero exit status")
    if want_help
        # TODO: we throw a error here, because most RNAstructure
        # programs return a non-zero exit status when help is
        # requested, but `partition` doesn't
        error("help string requested")
    end
    ct_structs = parse_ct_format(res)
    dbns = [pairtable_to_dbn(pt) for (_, _, pt) in ct_structs]
    return dbns
end

"""
    subopt(seq; [verbose, args]) -> subopt_en_structures

Calculate suboptimal structures for an RNA sequence `seq`. Returns a
vector of secondary structures in dot-bracket notation and their
energies.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for
details on command-line arguments that can be passed as `args`.

In particular, the `--maximum`, `--percent`, and `--window` arguments
can be used to increase or decrease the number of secondary structures
returned.
"""
function subopt(seq::AbstractString; verbose::Bool=false, args=``)
    if ("-h" in args) || ("--help" in args)
        exitcode, res, out, err = run_Fold(""; args=`-h`)
        println(out, err)
        error("help string requested")
    end
    exitcode, res, out, err = run_Fold(seq; args)
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
    subopt_all(seq; [verbose, args]) -> subopt_en_structures

Generate all suboptimal structures for an RNA sequence `seq`. Returns
a vector of secondary structures in dot-bracket notation and their
energies.

See the [RNAstructure AllSub
documentation](https://rna.urmc.rochester.edu/Text/AllSub.html) for
details on command-line arguments that can be passed as `args`.
"""
function subopt_all(seq::AbstractString; verbose::Bool=false, args=``)
    if ("-h" in args) || ("--help" in args)
        exitcode, res, out, err = run_AllSub(""; args=`-h`)
        println(out, err)
        error("help string requested")
    end
    exitcode, res, out, err = run_AllSub(seq; args)
    if verbose || exitcode != 0
        println("stdout of AllSub:")
        println(out)
        println("stderr of AllSub:")
        println(err)
    end
    # TODO: make energy parsing a global helper fn
    get_energy(str) = try parse(Float64, split(str, "=")[2]) catch; 0.0 end * u"kcal/mol"
    exitcode == 0 || error("AllSub returned non-zero exit status ($exitcode)")
    ct_structs = parse_ct_format(res)
    en_dbns = [(pairtable_to_dbn(pt), get_energy(title)) for (title, _, pt) in ct_structs]
    return en_dbns
end

"""
    run_AllSub(seq; [args]) -> exitcode, res, out, err

Run the `AllSub` program from RNAstructure.

See the [RNAstructure AllSub
documentation](https://rna.urmc.rochester.edu/Text/AllSub.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_AllSub(seq::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do seqpath, _
            _write_dbn_fasta(seqpath, seq)
            cmd = `$(RNAstructure_jll.AllSub()) $seqpath $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_ct2dot(ct::AbstractString, i::Integer=1; [args]) -> exitcode, res, out, err

Run the `ct2dot` program from RNAstructure. This converts structure
`i` from `ct` into dot-bracket notation.

See the [RNAstructure ct2dot
documentation](https://rna.urmc.rochester.edu/Text/ct2dot.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_ct2dot(ct::AbstractString, i::Integer=1; args=``)
    i >= 1 || throw(ArgumentError("structure number must be >= 1"))
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do ctpath, _
            write(ctpath, ct)
            cmd = `$(RNAstructure_jll.ct2dot()) $ctpath $i $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_CycleFold(seq; [args]) -> exitcode, out, err

Run the `CycleFold` program from RNAstructure.

See the [RNAstructure CycleFold
documentation](https://rna.urmc.rochester.edu/Text/CycleFold.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_CycleFold(seqs::Vector{<:AbstractString}; args=``)
    want_turbo = "-t" in args || "-T" in args || "--turbo" in args
    exitcode = 0
    res = out = err = ""
    mktemp() do seqpath, _
        if want_turbo
            # Note: CycleFold (RNAstructure-6.4.0) needs a title
            # string in the fasta file when using --turbo, otherwise
            # empty results and no error message are returned
            _write_fasta(seqpath, [("seq" * string(i), s) for (i,s) in enumerate(seqs)])
        else
            _write_fasta(seqpath, seqs)
        end
        cmd = `$(RNAstructure_jll.CycleFold()) $seqpath $args`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
end

run_CycleFold(seq::AbstractString; args=``) = run_CycleFold([seq]; args)

"""
    run_dot2ct(seq, dbn; [args]) -> exitcode, res, out, err
    run_dot2ct(dbn; [args]) -> exitcode, res, out, err

Run the `dot2ct` program from RNAstructure.

See the [RNAstructure dot2ct
documentation](https://rna.urmc.rochester.edu/Text/dot2ct.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_dot2ct(seq::AbstractString, dbn::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do dbnpath, _
            _write_dbn_fasta(dbnpath, seq, dbn)
            cmd = `$(RNAstructure_jll.dot2ct()) $dbnpath $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

run_dot2ct(dbn::AbstractString; args=``) =
    run_dot2ct("N"^length(dbn), dbn; args)

"""
    run_draw(dbn, [seq]; [args]) -> exitcode, res, out, err

Run the `draw` program from RNAstructure.

See the [RNAstructure draw
documentation](https://rna.urmc.rochester.edu/Text/draw.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_draw(dbn::AbstractString, seq::AbstractString;
                  args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do dbnpath, _
            _write_dbn_fasta(dbnpath, seq, dbn)
            cmd = `$(RNAstructure_jll.draw()) $dbnpath $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_EDcalculator(seq, dbn; [args]) -> exitcode, out, err
    run_EDcalculator(seq, dbns; [args])

Run the `EDcalculator` program from RNAstructure.

See the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html)
for details on command-line arguments that can be passed as `args`.
"""
function run_EDcalculator(seq::AbstractString, dbn::AbstractString;
                          args=``)
    return run_EDcalculator(seq, [dbn]; args)
end

function run_EDcalculator(seq::AbstractString, dbns::Vector{<:AbstractString};
                          args=``)
    exitcode = 0
    out = err = ""
    mktemp() do dbnpath, _
        _write_dbn_fasta(dbnpath, seq, dbns)
        cmd = `$(RNAstructure_jll.EDcalculator()) $dbnpath $args`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
end

"""
    run_efn2(seq, dbn; [args]) -> exitcode, res, out, err
    run_efn2(seq, dbns; [args])

Run the `efn2` program from RNAstructure. Returns the exitcode of the
`efn2` program, the contents of the results file `res`, the stdout
`out`, and stderr `err` output as Strings.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for
details on command-line arguments that can be passed as `args`.
"""
run_efn2(seq::AbstractString, dbn::AbstractString; args=``) =
    run_efn2(seq, [dbn]; args)

function run_efn2(seq::AbstractString,
                  dbns::Vector{<:AbstractString};
                  args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do dbnpath, _
        mktemp() do respath, _
            _write_dbn_fasta(dbnpath, seq, dbns)
            cmd = `$(RNAstructure_jll.efn2()) $dbnpath $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_EnsembleEnergy(seq; [args]) -> exitcode, out, err

Run the `EnsembleEnergy` program from RNAstructure. Returns the
exitcode of the `EnsembleEnergy` program, the stdout `out`, and stderr
`err` output as Strings.

See the [RNAstructure EnsembleEnergy
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html)
for details on command-line arguments that can be passed as `args`.
"""
function run_EnsembleEnergy(seq::AbstractString; args=``)
    exitcode = 0
    out = err = ""
    mktemp() do seqpath, _
        _write_dbn_fasta(seqpath, seq)
        cmd = `$(RNAstructure_jll.EnsembleEnergy()) $seqpath --sequence $args`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
end

"""
    run_Fold(seq; [args]) -> exitcode, res, out, err

Run the `Fold` program from RNAstructure.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_Fold(seq::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do seqpath, _
            _write_dbn_fasta(seqpath, seq)
            cmd = `$(RNAstructure_jll.Fold()) $seqpath $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_MaxExpect(seq; [args]) -> exitcode, res, out, err

Run the `MaxExpect` program from RNAstructure.

See the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html) for
details on command-line arguments that can be passed as `args`.
"""
function run_MaxExpect(pf_savefile::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        cmd = `$(RNAstructure_jll.MaxExpect()) $pf_savefile $respath $args`
        exitcode, out, err = _runcmd(cmd)
        res = read(respath, String)
    end
    return exitcode, res, out, err
end

"""
    run_partition!(pf_savefile, seq; [args]) -> exitcode, out, err

Run the `partition` program from RNAstructure. The partition function
save file is saved to the path given by `pf_savefile`.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html)
for details on command-line arguments that can be passed as `args`.
"""
function run_partition!(pf_savefile::AbstractString, seq::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do seqpath, _
        _write_dbn_fasta(seqpath, seq)
        cmd = `$(RNAstructure_jll.partition()) $seqpath $pf_savefile $args`
        exitcode, out, err = _runcmd(cmd)
    end
    return exitcode, out, err
end

"""
    run_ProbabilityPlot(pf_savefile; [args]) -> exitcode, res, out, err

Run the `ProbabilityPlot` program from RNAstructure.  The partition
function save file `pf_savefile` is the path to a file previously
generated with `run_partition!`.

See the [RNAstructure ProbabilityPlot
documentation](https://rna.urmc.rochester.edu/Text/ProbabilityPlot.html)
for details on command-line arguments that can be passed as `args`.
"""
function run_ProbabilityPlot(pf_savefile::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        cmd = `$(RNAstructure_jll.ProbabilityPlot()) $pf_savefile $respath $args`
        exitcode, out, err = _runcmd(cmd)
        res = read(respath, String)
    end
    return exitcode, res, out, err
end

"""
    run_RemovePseudoknots(seq, dbn; [args]) -> exitcode, res, out, err

Run the `RemovePseudoknots` program from RNAstructure.

See the [RNAstructure RemovePseudoknots
documentation](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html)
for details on command-line arguments that can be passed as `args`.
"""
function run_RemovePseudoknots(seq::AbstractString, dbn::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do ctpath, _
            ps, ct, _... = run_dot2ct(seq, dbn)
            ps == 0 || error("dot2ct returned non-zero exit status")
            open(ctpath, "w") do io
                write(io, ct)
            end
            cmd = `$(RNAstructure_jll.RemovePseudoknots()) $ctpath $respath $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

"""
    run_stochastic(seq; [args]) -> exitcode, res, out, err

Run the `stochastic` program from RNAstructure.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html)
for details on command-line arguments that can be passed as `args`.
"""
function run_stochastic(seq::AbstractString; args=``)
    exitcode = 0
    res = out = err = ""
    mktemp() do respath, _
        mktemp() do seqpath, _
            _write_dbn_fasta(seqpath, seq)
            cmd = `$(RNAstructure_jll.stochastic()) $seqpath $respath --sequence $args`
            exitcode, out, err = _runcmd(cmd)
            res = read(respath, String)
        end
    end
    return exitcode, res, out, err
end

end # module RNAstructure
