function partial_from_marginal_error(T, S)
    ArgumentError(
        """Cannot compute partial $(get_base_estimate_name(T)) from marginal $(get_base_estimate_name(S)), either call `k_function`
        or pass the original spectrum or data if you want the partial version."""
    )
end
