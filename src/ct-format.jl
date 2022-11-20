# CT (connectivity table, .ct) format for nucleic acid secondary structure
# - can describe pseudoknots, multiple strands, circular strands,
#   and record natural numbering (e.g. numbering from the
#   literature that deviates from consecutive numbering)
# - each site can interact with either one or no other sites (basepairs)
# - could in theory be used for proteins and other linear polymers
# Reference: https://rna.urmc.rochester.edu/Text/File_Formats.html#CT
# Format description:
# - header line: total_number_of_bases title_of_sequence
# - other lines:
#   - base_number(index n)
#   - base(A,C,G,T,U,X)
#   - prev_index(n-1)
#   - next_index(n+1)
#   - base_pairing_partner_index(0 if unpaired)
#   - natural numbering
# There can be multiple structures in one file, each preceded by
# its own header line.

# TODO
# This is the same function as in LinearFold.jl, put this in some
# BioRNAformats.jl package maybe?
#
# Note: This implementation has improvements that are not in the
# implementation in LinearFold.jl
# - doesn't fail if title is not exactly one word
# - returns title
# - returns sequence
# - can read multiple structures contained in one file

"""
    parse_ct_format(ctstr) -> [title, sequence, pairtable]

Parse secondary structures in ct (connectivity table, .ct) format.
"""
function parse_ct_format(ctstr::AbstractString)
    # TODO
    # - support multiple strands in one structure (iprev, inext)
    # - support circular strands
    iobuf = IOBuffer(ctstr)
    re_emptyline = r"^\s*$"
    re_header = r"^\s*(\d+)(\s*$|\s+(.*)$)"
    results = Tuple{String,Vector{String},Vector{Int}}[]
    want_header = true
    n = 0
    title = ""
    seq = String[]
    pt = Int[]
    for line in eachline(iobuf)
        # skip empty lines
        if occursin(re_emptyline, line)
            continue
        end
        if want_header
            # parse header line
            m = match(re_header, line)
            if isnothing(m)
                error("error parsing header line: $line")
            end
            n_str = m.captures[1]
            title = if isnothing(m.captures[3])
                ""
            else
                String(m.captures[3])
            end
            n = parse(Int, n_str)
            pt = zeros(Int, n)
            seq = Vector{String}(undef, n)
            want_header = false
            continue
        end
        # parse line with base-pairing information
        a = split(line)
        # we ignore extra entries at the end, as some programs put comments there
        if length(a) < 6
            error("Error: not enough entries (has to be >= 6) in line: $line")
        end
        i_str, base_str, iprev_str, inext_str, j_str, natural_idx_str = a[1:6]
        i = parse(Int, i_str)
        iprev = parse(Int, iprev_str)
        inext = parse(Int, inext_str)
        j = parse(Int, j_str)
        natural_idx = String(natural_idx_str)
        pt[i] = j
        seq[i] = String(base_str)
        if i == n
            push!(results, (title, seq, pt))
            want_header = true
        end
    end
    return results
end

"""
    print_ct_format(io, pairtable; title, seq)

Print secondary structure given as `pairtable` to `io` in ct
(connectivity table, .ct) format. Optionally, the `title` and sequence
`seq` can be set.

Pairtable format
- pairtable[i] == 0 means position i is unpaired
- pairtable[i] == j means position i is paired to position j
"""
function print_ct_format(io::IO, pairtable::Vector{<:Integer};
                         title::AbstractString="", seq::AbstractString="")
    # TODO
    # - support circular strands
    # - support multiple strands
    # - support alternative sequence numbering
    n = length(pairtable)
    if seq == ""
        seq = "N"^n
    end
    if length(seq) != n
        error("pairtable and seq must have same length")
    end
    println(io, "$n $title")
    for i = 1:n
        j = pairtable[i]
        println(io, "$i  $(seq[i])  $(i-1)  $(i+1)  $j  $i")
    end
end
