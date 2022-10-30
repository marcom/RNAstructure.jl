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

Parse secondary structures in .ct (connectivity table) format.
"""
function parse_ct_format(ctstr::AbstractString)
    # TODO
    # - preserve whitespace in title string
    # - support multiple strands (iprev, inext), circular strands
    # - seq: efficiency of string concat
    # - seq: String or Vector{Char} ?
    #
    # Notes
    # .ct file format (connectivity table)
    # - can describe pseudoknots, multiple strands, circular strands,
    #   and record natural numbering (e.g. numbering from the
    #   literature that doesn't follow consecutive numbering)
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
    iobuf = IOBuffer(ctstr)
    re_emptyline = r"^\s*$"
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
            a = split(line)
            if length(a) == 0
                error("Error in header line of ct file, expected at least one entry, ",
                      "got $(length(a)). Line was: $line")
            end
            # TODO: not quite right, whitespace in title not preserved
            # -> capture with regex
            n_str = a[1]
            title = if length(a) > 1
                join(a[2:end], " ")
            else
                ""
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
