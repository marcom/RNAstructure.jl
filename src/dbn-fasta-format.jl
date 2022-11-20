# "dbn-fasta" file format
# format:
# >title        (on one line)
# SEQUENCE      (on one line)
# STRUCTURE1    (on one line)
# STRUCTURE<N>...

function _write_dbn_fasta(path::AbstractString, seq::AbstractString, dbn::AbstractString;
                          title::AbstractString="")
    return _write_dbn_fasta(path, seq, [dbn]; title)
end

function _write_dbn_fasta(path::AbstractString, seq::AbstractString;
                          title::AbstractString="")
    return _write_dbn_fasta(path, seq, String[]; title)
end

function _write_dbn_fasta(path::AbstractString, seq::AbstractString,
                          dbns::Vector{<:AbstractString}; title="")
    open(path, "w") do io
        _write_dbn_fasta(io, seq, dbns; title)
    end
end

function _write_dbn_fasta(io::IO, seq::AbstractString,
                          dbns::Vector{<:AbstractString}; title="")
    println(io, ">", title)
    println(io, seq)
    for dbn in dbns
        println(io, dbn)
    end
end

# returns: title, seq, [dbns]
function _parse_dbn_fasta(dbnfasta::AbstractString)
    lines = readlines(IOBuffer(dbnfasta))
    if length(lines) < 3
        throw(ArgumentError("expected at least 3 input lines, only found $(length(lines))"))
    end
    title_arr = collect(lines[1])
    seq = lines[2]
    dbns = lines[3:end]
    if title_arr[1] != '>'
        throw(ArgumentError("expected title line to start with '>'"))
    end
    title = join(title_arr[2:end])
    return title, seq, dbns
end
