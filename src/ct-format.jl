using FoldRNA: Pairtable, isunpaired

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
#
# TODO
# - support circular strands
# - support multiple strands
# - support alternative sequence numbering
function print_ct_format(io::IO, dbn::AbstractString;
                         seq::AbstractString="", title::AbstractString="")
    n = length(dbn)
    if seq == ""
        seq = "N"^n
    end
    println(io, "$n $title")
    pt = Pairtable(dbn)
    for i = 1:n
        j = isunpaired(pt, i) ? 0 : pt.pairs[i]
        println(io, "$i  $(seq[i])  $(i-1)  $(i+1)  $j  $i")
    end
end
