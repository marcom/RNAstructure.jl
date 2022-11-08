struct Plot
    svg::String
end

function Base.showable(mime::Type{T}, p::Plot) where {T <: MIME}
    if mime === MIME"image/svg+xml"
        return true
    end
    return false
end

function Base.show(io::IO, mime::MIME"image/svg+xml", p::Plot)
    println(io, p.svg)
end

"""
    plot(dbn, seq; [verbose, args]) -> Plot

Plots a secondary structure `dbn` and nucleotide sequence `seq` in SVG
format. The plot will be shown directly in Jupyter or Pluto notebooks.

See the [RNAstructure draw
documentation](https://rna.urmc.rochester.edu/Text/draw.html) for
details on command-line arguments that can be passed as `args`.
"""
function plot(dbn::AbstractString, seq::AbstractString;
              verbose::Bool=false, args=``)
    args = `--svg -n 1 $args`
    ps, res, out, err = run_draw(dbn, seq; args)
    if verbose
        println("stdout of draw\n", out, "stderr of draw\n", err)
    end
    if ps != 0 || res == ""
        error("error running draw\n", out, "\n", err)
    end
    return Plot(res)
end
