getextrainformation(::AbstractEstimate) = () # if you need additional information, override this method
getbaseestimatename(est) = getbaseestimatename(typeof(est))
function getbaseestimatename(T::Type{<:AbstractEstimate}) # please override this if you want a different name, you should call `getestimatename` when using this as some types may overload that in certain cases
    typename = string(nameof(T))
    # Convert CamelCase to lowercase words with spaces
    lowercase(join(split(typename, r"(?=[A-Z])"), " "))
end
function getshortbaseestimatename(T::Type{<:AbstractEstimate})
    getbaseestimatename(T::Type{<:AbstractEstimate})
end
