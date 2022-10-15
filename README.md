# RNAstructure.jl

Unofficial Julia interface to the
[RNAstructure](https://rna.urmc.rochester.edu/RNAstructure.html)
program suite for RNA structure prediction and analysis.

## Installation

```
] add https://github.com/marcom/RNAstructure.jl
```

## Usage

### Free energy of folding

```julia
using RNAstructure

# returns energy and experimental uncertainty
energy("GGGAAACCC",
       "(((...)))")

# pseudoknot
energy("GGGAAAAGGGAAAACCCAAAACCC",
       "(((....[[[....)))....]]]")

# set temperature to 300 K
energy("GGGAAAAGGGAAAACCCAAAACCC",
       "(((....[[[....)))....]]]";
       cmdline_opts=["-T", "300"])

# multiple structures, returns array of results
energy("GGGAAACCC",
      ["(((...)))",
       "((.....))"])

# to see the help string with command-line options
energy("", ""; cmdline_opts=["-h"])
```

Note: the `energy` function calls the `efn2` program from the
RNAstructure_jll package and parses its output.
