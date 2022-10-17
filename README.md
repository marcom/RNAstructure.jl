# RNAstructure.jl

Unofficial Julia interface to the
[RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html)
program suite for RNA structure prediction and analysis.

## Installation

```julia
using Pkg
pkg"add https://github.com/marcom/RNAstructure.jl"
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
RNAstructure.edcalculator("GGGAAACCC", "(((...)))")
RNAstructure.edcalculator("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.edcalculator("GGGAAACCC", "(((...)))"; cmdargs=`-a dna`)
```

### efn2

The `efn2` program calculates the folding free energy of a sequence
and one or more secondary structures.

See the [RNAstructure efn2
documentation](https://rna.urmc.rochester.edu/Text/efn2.html) for more
details and for command-line arguments that can be passed via
`cmdargs`.

```julia
RNAstructure.efn2("GGGAAACCC", "(((...)))")
RNAstructure.efn2("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.efn2("GGGAAACCC", "(((...)))"; cmdargs=`-T 300`)
```

### Fold

The `Fold` program calculates minimum free energy (mfe) and suboptimal
structures.

See the [RNAstructure Fold
documentation](https://rna.urmc.rochester.edu/Text/Fold.html) for more
details and for command-line arguments that can be passed as
`cmdargs`.

```julia
RNAstructure.fold("GGGAAACCC")
RNAstructure.fold("GGGAAACCC"; cmdargs=`-mfe`)
```

