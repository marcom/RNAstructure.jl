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

function plot(dbn::AbstractString, seq::AbstractString; cmdargs=``)
    ps, res, out, err = run_draw(dbn, seq; cmdargs=`--svg -n 1 $cmdargs`)
    if ps != 0 || res == ""
        error("error running draw\n", out, "\n", err)
    end
    return Plot(res)
end
