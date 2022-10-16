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
       cmdargs=["-T", "300"])

# multiple structures, returns array of results
energy("GGGAAACCC",
      ["(((...)))",
       "((.....))"])

# to see the help string with command-line options
energy("", ""; cmdargs="-h")

# efn2 program call without output parsing
RNAstructure.efn2("GGGAAACCC", "(((...)))")
```

Note: the `energy` function calls the `efn2` program from the
RNAstructure_jll package and parses its output.

### Sequence design

```julia
target = "(((...)))"

# returns designed sequence and random seed used for design
design(target)

# set the random number seed used by the design process
seed = 42
design(target; cmdargs=["-s", seed])

# show additional possible command-line args for cmdargs
design(""; cmdargs="-h")
```
