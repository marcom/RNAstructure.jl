# RNAstructure.jl

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


### Minimum free energy (MFE) and structure

The `mfe` function calculates the minimum free energy and the
corresponding minimum free energy structure of an RNA
sequence. Internally, this function calls the `Fold` program from
RNAstructure.

Additional information on the `Fold` program and possible command-line
options that can be passed via `cmdargs` can be found at the
[RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html).

```julia
# returns mfe and structure
mfe("GGGAAACCC")              # -> (-1.2 kcal mol^-1, "(((...)))")

# set temperature to 300 K
mfe(seq; cmdargs=`-T 300`)    # -> (-1.9 kcal mol^-1, "(((...)))")

# show possible options for cmdargs
mfe(""; cmdargs=`-h`)
```


### Free energy of folding

The `energy` function calls the `efn2` program and parses its
output. It calculates the folding free energy and experimental
uncertainty of a sequence and one or more secondary structures.

Additional information on the `efn2` program and possible command-line
options that can be passed via `cmdargs` can be found at the
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
       cmdargs=`-T 300`)

# multiple structures, returns array of results
energy("GGGAAACCC",
      ["(((...)))",
       "((.....))"])

# show possible options for cmdargs
energy("", ""; cmdargs=`-h`)
```


### Basepair probabilities

The `bpp` function calls the `partition` and `ProbabilityPlot`
programs from RNAstructure to calculate the basepair probabilities for
an RNA sequence.

```julia
bpp("GGGAAACCC")  # -> 9x9 Matrix

# show possible options for cmdargs
bpp(""; cmdargs=`-h`)
```


### Sequence design

The `design` function calls the `design` program from RNAstructure.

Additional information on the `design` program and possible
command-line options that can be passed via `cmdargs` can be found at
the [RNAstructure design
documentation](https://rna.urmc.rochester.edu/Text/design.html).

```julia
target = "(((...)))"

# returns designed sequence and random seed used for design
design(target)

# set the random number seed used by the design process
seed = 42
design(target; cmdargs=`-s $seed`)

# show possible options for cmdargs
design(""; cmdargs=`-h`)
```


### Ensemble defect

The `ensemble_defect` function calls the `EDcalculator` program from
RNAstructure. It calculates the ensemble defect and normalised
ensemble defect of a sequence and one or more secondary structures.

Additional information on the `EDcalculator` program and possible
command-line options that can be passed via `cmdargs` can be found at
the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html).

```julia
seq = "GGGAAACCC"
dbn = "(((...)))"
dbns = [dbn, "((.....))"]
ensemble_defect(seq, dbn)
ensemble_defect(seq, dbns)
ensemble_defect("AAACCCTTT", "(((...)))"; cmdargs=`-a dna`)

# show possible options for cmdargs
ensemble_defect("", ""; cmdargs=`-h`)
```

### Convert dot-bracket notation to ct format

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

dbn2ct("GGGAAACCC", "(((...)))")
```

## Basic API to RNAstructure programs

These functions setup input files automatically and read output files,
but don't parse the results.  They typically return the exit status of
the RNAstructure program, the contents of the output file, and
stdout/stderr output. Additional command-line arguments can be passed
to the programs with the keyword argument `cmdargs`.

### EDcalculator

The `EDcalculator` program calculates the ensemble defect of a
sequence and one or more secondary structures.

See the [RNAstructure EDcalculator
documentation](https://rna.urmc.rochester.edu/Text/EDcalculator.html)
for more details and for command-line arguments that can be passed via
`cmdargs`.

```julia
RNAstructure.run_EDcalculator("GGGAAACCC", "(((...)))")
RNAstructure.run_EDcalculator("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.run_EDcalculator("GGGAAACCC", "(((...)))"; cmdargs=`-a dna`)
```

### efn2

The `efn2` program calculates the folding free energy of a sequence
and one or more secondary structures.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for more
details and for command-line arguments that can be passed via
`cmdargs`.

```julia
RNAstructure.run_efn2("GGGAAACCC", "(((...)))")
RNAstructure.run_efn2("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.run_efn2("GGGAAACCC", "(((...)))"; cmdargs=`-T 300`)
```

### EnsembleEnergy

The `EnsembleEnergy` program calculates the ensemble energy of
structures for an RNA sequence, given by the formula `-RT log(Q)`.

See the [RNAstructure EnsembleEnergy
documentation](https://rna.urmc.rochester.edu/Text/EnsembleEnergy.html)
for more details and for command-line arguments that can be passed via
`cmdargs`.

```julia
RNAstructure.run_EnsembleEnergy("GGGAAACCC")
```

### Fold

The `Fold` program calculates minimum free energy (mfe) and suboptimal
structures.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for more
details and for command-line arguments that can be passed as
`cmdargs`.

```julia
RNAstructure.run_Fold("GGGAAACCC")
RNAstructure.run_Fold("GGGAAACCC"; cmdargs=`-mfe`)
```

### MaxExpect

The `MaxExpect` program predicts the maximum expected accuracy (MEA)
structure for an RNA sequence.

See the [RNAstructure MaxExpect
documentation](https://rna.urmc.rochester.edu/Text/MaxExpect.html) for
more details and for command-line arguments that can be passed via
`cmdargs`.

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
`cmdargs`.

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
`cmdargs`.

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
`cmdargs`.

```julia
# maximise basepairs in returned structure
dbn = "((...[[[[...))..]].]]"
RNAstructure.run_RemovePseudoknots("N"^length(dbn), dbn; cmdargs=`-m`)

# return pseudoknot-free structure with lowest folding free energy at
# a temperature of 300 K for a given sequence
seq = "GGAAAAUGCAAACCAAGCAAU"
RNAstructure.run_RemovePseudoknots(seq, dbn; cmdargs=`-T 300`)
```

### stochastic

The `stochastic` program samples from the Boltzmann ensemble of
secondary structures.

See the [RNAstructure stochastic
documentation](https://rna.urmc.rochester.edu/Text/stochastic.html)
for more details and for command-line arguments that can be passed via
`cmdargs`.

```julia
RNAstructure.run_stochastic("GGGAAACCC")
```

## Related Julia packages

- [ViennaRNA.jl](https://github.com/marcom/ViennaRNA.jl)
- [LinearFold.jl](https://github.com/marcom/LinearFold.jl)
- [PlotRNA.jl](https://github.com/marcom/PlotRNA.jl)
