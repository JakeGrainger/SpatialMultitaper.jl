# some types for determining when to compute split partial spectra

abstract type PartialType end
struct UsualPartial <: PartialType end
struct SplitPartial <: PartialType end