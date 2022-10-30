# TODO
# This is the same function as in LinearFold.jl, put this in some
# BioRNAformats.jl package maybe?
#
# Note: This implementation has improvements that are not in the
# implementation in LinearFold.jl
# - doesn't fail if title is not exactly one word
# - returns title
# - returns sequence

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
    # - first line: total_number_of_bases title_of_sequence
    # - other lines:
    #   - base_number(index n)
    #   - base(A,C,G,T,U,X)
    #   - prev_index(n-1)
    #   - next_index(n+1)
    #   - base_pairing_partner_index(0 if unpaired)
    #   - natural numbering
    iobuf = IOBuffer(ctstr)
    firstline = readline(iobuf)
    a = split(firstline)
    if length(a) == 0
        error("Error in first line of ct file, expected at least one entry, got $(length(a)). Line was: $firstline")
    end
    n_str = a[1]
    title = if length(a) > 1
        join(a[2:end], " ")  # TODO: not quite right, whitespace not preserved
    else
        ""
    end
    n = parse(Int, n_str)
    seq = ""
    pt = zeros(Int, n) # pairtable
    for line in eachline(iobuf)
        if isempty(strip(line))
            continue
        end
        a = split(line)
        # we ignore extra entries at the end, as some programs put comments there
        if length(a) < 6
            error("Error: not enough entries (has to be >= 6) in line: $line")
        end
        i_str, base_str, iprev_str, inext_str, j_str, natural_idx_str = a[1:6]
        i = parse(Int, i_str)
        seq *= String(base_str)
        iprev = parse(Int, iprev_str)
        inext = parse(Int, inext_str)
        j = parse(Int, j_str)
        natural_idx = String(natural_idx_str)
        pt[i] = j
    end
    return title, seq, pt
end
