# RNAstructure.jl

Unofficial Julia interface to the RNAstructure program suite for RNA
structure prediction and analysis.

## Installation

```
] add https://github.com/marcom/FoldRNA.jl
] add https://github.com/marcom/RNAstructure.jl
```

## Usage

```julia
using RNAstructure

energy("GGGAAACCC",
       "(((...)))")

# pseudoknot
energy("GGGAAAAGGGAAAACCCAAAACCC",
       "(((....[[[....)))....]]]")

# set temperature to 300 K
energy("GGGAAAAGGGAAAACCCAAAACCC",
       "(((....[[[....)))....]]]";
       cmdline_opts=["-T", "300"])

# to see the help string with command-line options
energy("A", "."; cmdline_opts=["-h"])
```