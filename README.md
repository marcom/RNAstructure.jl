# RNAstructure.jl

[![Build Status](https://github.com/marcom/RNAstructure.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/marcom/RNAstructure.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Aqua QA](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

Unofficial Julia interface to the
[RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html)
program suite for RNA structure prediction and analysis.  Please cite
the appropriate publications listed on the RNAstructure website if you
use this library.

## Installation

Enter the package mode from the Julia REPL by pressing `]` and then
install with

```
add RNAstructure
```

## Usage

```julia
using RNAstructure
```

### Note: sequence conventions

Sequences passed to RNAstructure use the following convention:
- uppercase character: normal nucleotide, U equivalent to T
- lowercase character: nucleotide cannot form basepairs
- X or N character: unknown base or base that cannot interact with
  others (cannot pair or stack)

See the [RNAstructure manual section for
sequences](https://rna.urmc.rochester.edu/Text/File_Formats.html#Sequence)
for more details.

Some programs make exceptions to these rules, check the manual pages
of the RNAstructure programs for details on any differences.

### Note: Overriding energy parameter directories

The environment variables `RNASTRUCTURE_JL_DATAPATH` can be set to
override the directory where energy parameters are read from. For the
`cyclefold_*` functions the environment variable is called
`RNASTRUCTURE_JL_CYCLEFOLD_DATAPATH`.

In the original RNAstructure program these environment variables are
called `DATAPATH` and `CYCLEFOLD_DATAPATH`. `RNAstructure.jl` (this
package) sets these environment variables automatically to the
corresponding installation directory of the `RNAstructure_jll` binary
package.  The names of the env vars were changed to avoid clashes with
possible settings you might already have in your shell startup files
from a pre-existing manual RNAstructure installation, which could be a
different version and have different parameters. In this way, you can
be sure that this package uses the correct parameters, while still
allowing to override them if necessary.

### Minimum free energy (MFE) and structure

The `mfe` function calculates the minimum free energy and the
corresponding minimum free energy structure of an RNA
sequence. Internally, this function calls the `Fold` program from
RNAstructure.

Additional information on the `Fold` program and possible command-line
options that can be passed via `args` can be found at the
[RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html).

```julia
# returns mfe and structure
mfe("GGGAAACCC")              # -> (-1.2 kcal mol^-1, "(((...)))")

# set temperature to 300 K
mfe(seq; args=`-T 300`)    # -> (-1.9 kcal mol^-1, "(((...)))")

# show possible options for args
mfe(""; args=`-h`)
```


### Suboptimal structures

Generate suboptimal structures for a nucleic acid
sequence. Internally, this function calls the `Fold` program from
RNAstructure.

Additional information on the `Fold` program and possible command-line
options that can be passed via `args` can be found at the
[RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html).

```julia
subopt("GGGAAACCC")
subopt("GGGGAAACCCC"; args=`-w 0 -p 100`)

# show possible options for args
subopt(""; args=`-h`)
```


### All suboptimal structures in an energy range

Generate all suboptimal structures in an energy range for a nucleic
acid sequence using the `AllSub` program from RNAstructure.

Additional information on the `AllSub` program and possible
command-line options that can be passed via `args` can be found at
the [RNAstructure AllSub
documentation](https://rna.urmc.rochester.edu/Text/AllSub.html).

```julia
subopt_all("GGGAAACCC")

# maximum absolute energy difference of 10 kcal/mol to the MFE, up to
# 500 percent relative difference to MFE
subopt_all("GGGGAAACCCC"; args=`-a 10 -p 500`)

# set temperature to 300 K
subopt_all("GGGGAAACCCC"; args=`-T 300`)

# show possible options for args
subopt_all(""; args=`-h`)
```


### Partition function (ensemble energy)

The `partfn` function calculates the partition function and returns
the ensemble free energy for a nucleotide sequence.

Additional information on the `EnsembleEnergy` program and possible
command-line options that can be passed via `args` can be found at
the [RNAstructure EnsembleEnergy
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html).

```julia
partfn("GGGAAACCC")

partfn("GGGAAACCC"; args=`--DNA`)

# show possible options for args_partition, args_maxexpect
partfn(""; args=`-h`)
```


### Probability of a structure

The `prob_of_structure` function calculates the probability of a
secondary structure for a given nucleotide sequence.

The supported args are those common to `energy` and `partfn`.

```julia
prob_of_structure("GGGAAACCC", "(((...)))")
```


### Maximum expected accuracy (MEA) structure

The `mea` function predicts the maximum expected accuracy structure
(and possibly suboptimals) for a nucleotide sequence.

Additional information on the `partition` program and possible
command-line options that can be passed via `args_partition` can be
found at the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html).

Additional information on the `MaxExpect` program and possible
command-line options that can be passed via `args_maxexpect` can be
found at the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html).

```julia
mea("GGGAAACCC")

mea("GGGAAACCC"; args_partition=`-T 300`, args_maxexpect=`-s 10 -w 0`)

# show possible options for args_partition, args_maxexpect
mea(""; args_partition=`-h`)
```


### Free energy of folding

The `energy` function calls the `efn2` program and parses its
output. It calculates the folding free energy and experimental
uncertainty of a sequence and one or more secondary structures.

Additional information on the `efn2` program and possible command-line
options that can be passed via `args` can be found at the
[RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html).

```julia
# returns energy and experimental uncertainty
energy("GGGAAACCC",
       "(((...)))")

# pseudoknot
energy("GGGAAAAGGGAAAACCCAAAACCC",
       "(((....[[[....)))....]]]")

# set temperature to 300 K
energy("GGGAAAAGGGAAAACCCAAAACCC",
       "(((....[[[....)))....]]]";
       args=`-T 300`)

# multiple structures, returns array of results
energy("GGGAAACCC",
      ["(((...)))",
       "((.....))"])

# show possible options for args
energy("", ""; args=`-h`)
```


### Basepair probabilities

The `bpp` function calls the `partition` and `ProbabilityPlot`
programs from RNAstructure to calculate the basepair probabilities for
an RNA sequence.

```julia
bpp("GGGAAACCC")  # -> 9x9 Matrix

# show possible options for args
bpp(""; args=`-h`)
```


### Sampling structures

Sample secondary structures from the Boltzmann ensemble of secondary
structures.

Additional information on the `stochastic` program and possible
command-line options that can be passed via `args` can be found at
the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html).

```julia
# returns a 1000-element Vector{String}
sample_structures("GGGAAACCC")

# show possible options for args
sample_structures(""; args=`-h`)
```


### Nucleotide cyclic motif model (CycleFold)

The `cyclefold_*` functions call the `CycleFold` program from
RNAstructure, which uses the nucleotide cyclic motif model by
(Parisien & Major, 2008).  This model allows for non-canonical and
canonical basepairs.

NOTE: use the energy with caution --- i think the energy unit is
kJ/mol, but i am not sure.

Additional information on the `CycleFold` program and possible
command-line options that can be passed via `args` can be found at the
[RNAstructure CycleFold
documentation](https://rna.urmc.rochester.edu/Text/CycleFold.html).

```julia
cyclefold_mea("GGGAAACCC")  # -> [9, 8, 7, 6, 0, 4, 3, 2, 1]
cyclefold_mfe("GGGAAACCC")  # -> (-7.8305 kJ mol^-1, [9, 8, 7, 6, 0, 4, 3, 2, 1])
cyclefold_bpp("GGGAAACCC")  # -> 9×9 Matrix{Float64}

# show possible options for args
cyclefold_mea(""; args=`-h`)
```


### Sequence design

The `design` function calls the `design` program from RNAstructure.

Additional information on the `design` program and possible
command-line options that can be passed via `args` can be found at
the [RNAstructure design
documentation](https://rna.urmc.rochester.edu/Text/design.html).

```julia
target = "(((...)))"

# returns designed sequence and random seed used for design
design(target)

# set the random number seed used by the design process
seed = 42
design(target; args=`-s $seed`)

# show possible options for args
design(""; args=`-h`)
```


### Ensemble defect

The `ensemble_defect` function calls the `EDcalculator` program from
RNAstructure. It calculates the ensemble defect and normalised
ensemble defect of a sequence and one or more secondary structures.

Additional information on the `EDcalculator` program and possible
command-line options that can be passed via `args` can be found at
the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html).

```julia
seq = "GGGAAACCC"
dbn = "(((...)))"
dbns = [dbn, "((.....))"]
ensemble_defect(seq, dbn)
ensemble_defect(seq, dbns)
ensemble_defect("AAACCCTTT", "(((...)))"; args=`-a dna`)

# show possible options for args
ensemble_defect("", ""; args=`-h`)
```


### Remove pseudoknots

The `remove_pseudoknots` function returns the pseudoknot-free
substructure with the maximum possible basepairs.

```julia
remove_pknots("(((...[[[[...)))...]]]]")  # -> "......((((.........))))"
```


### dbn2ct: convert dot-bracket notation to ct format

This function uses the `dot2ct` program from RNAstructure to convert a
secondary structure in dot-bracket notation and optionally a sequence
to the ct (connectivity table) format.

```julia
# if no sequence is given, it will be all 'N' in the resulting ct
# format output
dbn2ct("(((...)))")

# pseudoknots work as well
dbn2ct("(((...[[[...)))...]]]")
dbn2ct("(((...[[[...{{{...<<<...)))...]]]...}}}...>>>")

dbn2ct("(((...)))"; seq="GGGAAACCC")
dbn2ct(["(((...)))", "........."]; title="A sequence", seq="GGGAAACCC")
```

### ct2dbn: convert ct format to dot-bracket notation

This function uses the `ct2dot` program from RNAstructure to convert a
secondary structure and sequence in ct (connectivity table) format to
dot-bracket notation.

```julia
ct = dbn2ct("(((...)))"; title="TITLE", seq="GGGAAACCC")
print(ct)
ct2dbn(ct)  # -> (title = "TITLE", seq = "GGGAAACCC", dbn = "(((...)))")

ct2 = dbn2ct("........."; title="TITLE2", seq="NNNAAANNN")
ct_twostruct = ct * ct2
print(ct_twostruct)
ct2dbn(ct_twostruct, 2)  # -> (title = "TITLE2", seq = "NNNAAANNN", dbn = ".........")
```


### Plotting a secondary structure

This function uses the `draw` program from RNAstructure to plot a
secondary structure in dot-bracket notation to SVG format.  This
should show an image when used in Jupyter and Pluto notebooks.

Additional information on the `draw` program and possible command-line
options that can be passed via `args` can be found at the
[RNAstructure draw
documentation](https://rna.urmc.rochester.edu/Text/draw.html).

```julia
plot("(((...)))", "GGGAAACCC")
plot("(((...)))", "GGGAAACCC"; args=`--circle`)
plot("(((...)))", "GGGAAACCC"; args=`--flat`)
plot("(((...)))", "GGGAAACCC"; args=`--uncircled`)
```


## Basic API to RNAstructure programs

These functions setup input files automatically and read output files,
but don't parse the results.  They typically return the exit status of
the RNAstructure program, the contents of the output file, and
stdout/stderr output. Additional command-line arguments can be passed
to the programs with the keyword argument `args`.

### AllSub

The `AllSub` program calculates all suboptimal structures within a
certain energy range.

See the [RNAstructure AllSub
documentation](https://rna.urmc.rochester.edu/Text/AllSub.html) for
more details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_AllSub("GGGAAACCC")
RNAstructure.run_AllSub("GGGAAACCC"; args=`-a 10 -p 500`)
```

### ct2dot

The `ct2dot` converts secondary structures in connectivity table (ct)
format to dot-bracket notation.

See the [RNAstructure ct2dot
documentation](https://rna.urmc.rochester.edu/Text/ct2dot.html) for
more details and for command-line arguments that can be passed via
`args`.

```julia
ct = dbn2ct("(((...)))")
RNAstructure.run_ct2dot(ct)

ct2 = ct * dbn2ct(".........")
RNAstructure.run_ct2dot(ct2, 2)
```

### dot2ct

The `dot2ct` converts secondary structures in dot-bracket notation to
connectivity table (ct) format.

See the [RNAstructure dot2ct
documentation](https://rna.urmc.rochester.edu/Text/dot2ct.html) for
more details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_dot2ct("(((...)))")
RNAstructure.run_dot2ct("(((...)))"; seq="GGGAAACCC")
RNAstructure.run_dot2ct(["(((...)))", "........."];
                        title="A sequence", seq="GGGAAACCC")
```

### draw

The `draw` program draws secondary structure diagrams.

See the [RNAstructure draw
documentation](https://rna.urmc.rochester.edu/Text/draw.html) for more
details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_draw("(((...)))", "GGGAAACCC"; args=`--svg`)
```

### EDcalculator

The `EDcalculator` program calculates the ensemble defect of a
sequence and one or more secondary structures.

See the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html)
for more details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_EDcalculator("GGGAAACCC", "(((...)))")
RNAstructure.run_EDcalculator("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.run_EDcalculator("GGGAAACCC", "(((...)))"; args=`-a dna`)
```

### efn2

The `efn2` program calculates the folding free energy of a sequence
and one or more secondary structures.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for more
details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_efn2("GGGAAACCC", "(((...)))")
RNAstructure.run_efn2("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.run_efn2("GGGAAACCC", "(((...)))"; args=`-T 300`)
```

### EnsembleEnergy

The `EnsembleEnergy` program calculates the ensemble energy of
structures for an RNA sequence, given by the formula `-RT log(Q)`.

See the [RNAstructure EnsembleEnergy
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html)
for more details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_EnsembleEnergy("GGGAAACCC")
```

### Fold

The `Fold` program calculates minimum free energy (mfe) and suboptimal
structures.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for more
details and for command-line arguments that can be passed as
`args`.

```julia
RNAstructure.run_Fold("GGGAAACCC")
RNAstructure.run_Fold("GGGAAACCC"; args=`-mfe`)
```

### MaxExpect

The `MaxExpect` program predicts the maximum expected accuracy (MEA)
structure for an RNA sequence.

See the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html) for
more details and for command-line arguments that can be passed via
`args`.

```julia
pf_savefile = "out.pfs"
RNAstructure.run_partition!(pf_savefile, "GGGAAACCC")
RNAstructure.run_MaxExpect(pf_savefile)
```

### partition

The `partition` program calculates the partition function and basepair
probabilities for an RNA sequence and saves this information in a
partition save file, which can then be used by other programs.

See the [RNAstructure partition
documentation](https://rna.urmc.rochester.edu/Text/partition.html) for
more details and for command-line arguments that can be passed via
`args`.

```julia
# write the partition function save file to "save.pfs", overwriting
# any data if the file already exists
RNAstructure.run_partition!("save.pfs", "GGGAAACCC")
```

### ProbabilityPlot

The `ProbabilityPlot` program extracts basepair probabilities from a
partition function save file generated with `partition` and can output
them as a text file or as a dot plot.

See the [RNAstructure ProbabilityPlot
documentation](https://rna.urmc.rochester.edu/Text/ProbabilityPlot.html)
for more details and for command-line arguments that can be passed via
`args`.

```julia
pf_savefile = "save.pfs"
RNAstructure.run_partition!(pf_savefile, "GGGAAACCC")
RNAstructure.run_ProbabilityPlot(pf_savefile)
```

### RemovePseudoknots

The `RemovePseudoknots` program removes pseudoknots from an RNA
secondary structure, returning either the structure with the most base
pairs or the structure with lowest folding free energy.

See the [RNAstructure RemovePseudoknots
documentation](https://rna.urmc.rochester.edu/Text/RemovePseudoknots.html)
for more details and for command-line arguments that can be passed via
`args`.

```julia
# maximise basepairs in returned structure
dbn = "((...[[[[...))..]].]]"
RNAstructure.run_RemovePseudoknots("N"^length(dbn), dbn; args=`-m`)

# return pseudoknot-free structure with lowest folding free energy at
# a temperature of 300 K for a given sequence
seq = "GGAAAAUGCAAACCAAGCAAU"
RNAstructure.run_RemovePseudoknots(seq, dbn; args=`-T 300`)
```

### stochastic

The `stochastic` program samples from the Boltzmann ensemble of
secondary structures.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html)
for more details and for command-line arguments that can be passed via
`args`.

```julia
RNAstructure.run_stochastic("GGGAAACCC")
```

## Related Julia packages

- [ViennaRNA.jl](https://github.com/marcom/ViennaRNA.jl)
- [LinearFold.jl](https://github.com/marcom/LinearFold.jl)
- [PlotRNA.jl](https://github.com/marcom/PlotRNA.jl)
