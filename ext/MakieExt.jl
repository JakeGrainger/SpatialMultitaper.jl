module MakieExt

import SpatialMultitaper: AbstractEstimate

using Makie

function Makie.convert_arguments(::Type{<:AbstractPlot}, est::AbstractEstimate)
    return collect(est)
end

end
