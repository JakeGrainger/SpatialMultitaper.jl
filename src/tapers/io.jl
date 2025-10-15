# Display and I/O methods for taper types.

function Base.show(io::IO, ::MIME"text/plain", taper::GridTaper)
    print(io, "GridTaper on $(size(taper.grid)) grid")
end

function Base.show(io::IO, ::MIME"text/plain", taper::InterpolatedTaper)
    print(
        io, "InterpolatedTaper (continuous, derived from $(size(taper.source_grid)) grid)")
end

function Base.show(io::IO, ::MIME"text/plain", taper::SinTaper{D}) where {D}
    modes_str = join(taper.modes, "×")
    print(io, "SinTaper{$D}(modes=$modes_str)")
end

# Abstract type fallbacks
function Base.show(io::IO, ::MIME"text/plain", taper::ContinuousTaper)
    print(io, "A continuous taper function.")
end

function Base.show(io::IO, ::MIME"text/plain", taper::DiscreteTaper)
    print(io, "A discrete taper function.")
end

function Base.show(io::IO, ::MIME"text/plain", family::TaperFamily)
    region_str = if family.concentration_region !== nothing
        if family.concentration_region isa Ball
            " (concentration region: Ball, radius=$(family.concentration_region.radius))"
        elseif family.concentration_region isa Box
            " (concentration region: Box)"
        else
            " (concentration region: $(typeof(family.concentration_region)))"
        end
    else
        ""
    end
    print(io, "A family of $(length(family)) taper functions$region_str")
end
