function Base.show(io::IO, ::MIME"text/plain", estimate::AbstractEstimate)
    D = embeddim(estimate)
    proc_names_1, proc_names_2 = processnames(estimate)
    println(io, getestimatename(estimate), " of a ", D, " dimensional process")
    println(io, "  between processes: ", proc_names_1, " and ", proc_names_2)
    println(io, "  evaluated at ", _printargument(estimate))
end

function _printargument(estimate::AbstractEstimate)
    _printargument(getargument(estimate))
end
function _printargument(arg::Tuple)
    join(arg, ", ")
end
_printargument(arg) = string(arg)
