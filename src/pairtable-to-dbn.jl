# TODO
# This function is similar to String(::Pairtable) in FoldRNA.jl, merge
# them later.

"""
    pairtable_to_dbn(pairtable::Vector{Int}) -> dbn::String

Convert a secondary structure in pairtable representation to a String
in dot-bracket notation.
"""
function pairtable_to_dbn(pairtable::Vector{Int})
    # TODO: handle pseudoknotted structures
    n = length(pairtable)
    dbn = Vector{Char}(undef, n)
    for i = 1:n
        j = pairtable[i]
        if j == 0
            dbn[i] = '.'
        else
            if j > i
                dbn[i] = '('
            else
                dbn[i] = ')'
            end
        end
    end
    return join(dbn)
end
