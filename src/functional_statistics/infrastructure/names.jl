"""
    get_base_estimate_name(est::AbstractEstimate) -> String
    get_base_estimate_name(T::Type{<:AbstractEstimate}) -> String

Get the base name of an estimate type for display purposes.

# Arguments
- `est::AbstractEstimate`: The estimate object
- `T::Type{<:AbstractEstimate}`: The estimate type

# Returns
- `String`: Human-readable name derived from the type name
"""
get_base_estimate_name(est) = get_base_estimate_name(typeof(est))

# Please override this if you want a different name. You should call `get_estimate_name`
# when using this as some types may overload that in certain cases.
function get_base_estimate_name(T::Type{<:AbstractEstimate})
    typename = string(nameof(T))
    # Convert CamelCase to lowercase words with spaces
    lowercase(join(split(typename, r"(?=[A-Z])"), " "))
end

"""
    get_short_base_estimate_name(T::Type{<:AbstractEstimate}) -> String
    get_short_base_estimate_name(est::AbstractEstimate) -> String

Get a shortened version of the base estimate name. Default implementation
returns the same as `get_base_estimate_name`.

# Arguments
- `T::Type{<:AbstractEstimate}`: The estimate type
- `est::AbstractEstimate`: The estimate object

# Returns
- `String`: Shortened human-readable name, primarily for conversion to R
"""
function get_short_base_estimate_name(T::Type{<:AbstractEstimate})
    get_base_estimate_name(T)
end
function get_short_base_estimate_name(est::AbstractEstimate)
    get_short_base_estimate_name(typeof(est))
end

"""
    get_estimate_name(T::Type{<:AbstractEstimate}) -> String
    get_estimate_name(est::AbstractEstimate) -> String

Get the full estimate name, including trait-specific prefixes (e.g., "partial").

# Arguments
- `T::Type{<:AbstractEstimate}`: The estimate type
- `est::AbstractEstimate`: The estimate object

# Returns
- `String`: Full estimate name with appropriate prefixes
"""
get_estimate_name(T::Type{<:MarginalAbstractEstimate}) = get_base_estimate_name(T)
function get_estimate_name(T::Type{<:PartialAbstractEstimate})
    "partial " * get_base_estimate_name(T)
end
get_estimate_name(est::AbstractEstimate) = get_estimate_name(typeof(est))

"""
    get_short_estimate_name(T::Type{<:AbstractEstimate}) -> String
    get_short_estimate_name(est::AbstractEstimate) -> String

Get the short estimate name, including trait-specific prefixes (e.g., "partial").

# Arguments
- `T::Type{<:AbstractEstimate}`: The estimate type
- `est::AbstractEstimate`: The estimate object

# Returns
- `String`: Short estimate name with appropriate prefixes
"""
get_short_estimate_name(T::Type{<:MarginalAbstractEstimate}) = get_short_base_estimate_name(T)
function get_short_estimate_name(T::Type{<:PartialAbstractEstimate})
    "partial " * get_short_base_estimate_name(T)
end
get_short_estimate_name(est::AbstractEstimate) = get_short_estimate_name(typeof(est))
