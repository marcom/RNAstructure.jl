# RNAstructure.jl

Unofficial Julia interface to the
[RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html)
program suite for RNA structure prediction and analysis.

## Installation

```
] add https://github.com/marcom/RNAstructure.jl
```

## Usage

```julia
using RNAstructure
```

Note: sequences passed to RNAstructure use the following conventions
- uppercase character: normal nucleotide, U equivalent to T
- lowercase character: nucleotide cannot form basepairs
- X or N character: unknown base or base that cannot interact with
  others (cannot pair or stack)

See the [RNAstructure manual section for
sequences](https://rna.urmc.rochester.edu/Text/File_Formats.html#Sequence)
for more details.


### Free energy of folding

The `energy` function calls the `efn2` program and parses its output.

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

# to see the help string with command-line options
energy("", ""; cmdargs=`-h`)
```


### Sequence design

```julia
target = "(((...)))"

# returns designed sequence and random seed used for design
design(target)

# set the random number seed used by the design process
seed = 42
design(target; cmdargs=`-s $seed`)

# show additional possible command-line args for cmdargs
design(""; cmdargs=`-h`)
```


## Basic API to RNAstructure programs

These functions setup input files automatically and read output files,
but don't parse the results.  They typically return the contents of
the output file, and stdout/stderr output.

### Fold

The `Fold` program calculates minimum free energy (mfe) and suboptimal
structures.

See the [RNAstructure Fold
manpage](https://rna.urmc.rochester.edu/Text/Fold.html) for details on
command-line arguments that can be passed as `cmdargs`.

```julia
RNAstructure.fold("GGGAAACCC")
RNAstructure.fold("GGGAAACCC"; cmdargs=`-mfe`)
```

### efn2

The `efn2` program calculates the folding free energy of a sequence
and one or more secondary structures.

See the [RNAstructure efn2
manpage](https://rna.urmc.rochester.edu/Text/efn2.html) for details on
command-line arguments that can be passed as `cmdargs`.

```julia
RNAstructure.efn2("GGGAAACCC", "(((...)))")
RNAstructure.efn2("GGGAAACCC", ["(((...)))", "((.....))"])
RNAstructure.efn2("GGGAAACCC", "(((...)))"; cmdargs=`-T 300`)
```
