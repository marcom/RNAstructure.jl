# fasta file format
# >title1
# SEQUENCE1
# >title2
# SEQUENCE2
# ...

function _write_fasta(path::AbstractString, title_seq_iterable)
    open(path, "w") do io
        _write_fasta(io, title_seq_iterable)
    end
end

function _write_fasta(io::IO, title_seq_iterable)
    for (title, seq) in title_seq_iterable
        println(io, ">", title, "\n", seq)
    end
end

_write_fasta(path::AbstractString, seqs::Vector{<:AbstractString}) =
    _write_fasta(path, [("", s) for s in seqs])

_write_fasta(path::AbstractString, seq::AbstractString) =
    _write_fasta(path, [("", seq)])

_write_fasta(io::IO, seqs::Vector{<:AbstractString}) =
    _write_fasta(io, [("", s) for s in seqs])

_write_fasta(io::IO, seq::AbstractString) =
    _write_fasta(io, [("", seq)])
