function partial_from_marginal_error(T, S)
    ArgumentError(
        """Cannot compute partial $(getbaseestimatename(T)) from marginal $(getbaseestimatename(S)), either call `k_function`
        or pass the original spectrum or data if you want the partial version."""
    )
end
