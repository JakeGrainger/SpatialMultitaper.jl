# Display and I/O methods for taper types.

function Base.show(io::IO, ::MIME"text/plain", taper::ContinuousTaper)
    print(io, "A continuous taper function.")
end

function Base.show(io::IO, ::MIME"text/plain", taper::DiscreteTaper)
    print(io, "A discrete taper function.")
end

function Base.show(io::IO, ::MIME"text/plain", taper::InterpolatedTaper)
    print(io, "A continuous interpolated taper function.")
end

function Base.show(io::IO, ::MIME"text/plain", taper::TaperFamily)
    print(io, "A family of $(length(taper)) taper functions.")
end
