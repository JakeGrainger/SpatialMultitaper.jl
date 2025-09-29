module TestUtils

"""
Main test utilities module for SpatialMultitaper.jl tests.

This module re-exports all test utilities from specialized submodules:
- TestData: Functions for generating test datasets
- NumericUtils: Numerical utilities for testing (reference implementations, etc.)
- TestHelpers: Additional helper functions for tests

Usage in test files:
```julia
include("test_utilities/TestUtils.jl")
using .TestUtils
```

This provides access to all test utilities in a clean, organized way.
"""

# Include all submodules
include("TestData.jl")
include("NumericUtils.jl")
include("TestHelpers.jl")

# Re-export everything from submodules
using .TestData
using .NumericUtils
using .TestHelpers

# Re-export all functions for convenient access
export make_points_example, make_grids_example, make_marked_example, make_mixed_example,
       basic_test_grid,  # from TestData
       slow_dft                                                                                            # from NumericUtils
# Add TestHelpers exports here when they exist

end # module TestUtils
