# Compact show method for use in arrays and other contexts
function Base.show(io::IO, estimate::AbstractEstimate)
    D = embeddim(estimate)
    proc_names_1, proc_names_2 = processnames(estimate)
    print(io, get_estimate_name(estimate), "(", D, "D, ", proc_names_1,
        "↔", proc_names_2, ", ", _printargument(estimate), ")")
end

# Detailed show method for standalone display
function Base.show(io::IO, ::MIME"text/plain", estimate::AbstractEstimate)
    D = embeddim(estimate)
    proc_names_1, proc_names_2 = processnames(estimate)
    println(io, get_estimate_name(estimate), " of a ", D, " dimensional process")
    println(io, "  between processes: ", proc_names_1, " and ", proc_names_2)
    println(io, "  evaluated at ", _printargument(estimate))
    println(io, "  with values of type ", typeof(get_estimates(estimate)))
end

function _printargument(estimate::AbstractEstimate)
    _printargument(get_evaluation_points(estimate))
end
function _printargument(arg::Tuple)
    join(arg, ", ")
end
_printargument(arg) = string(arg)
