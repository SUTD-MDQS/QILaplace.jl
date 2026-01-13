# Preamble file that contains the important functions and package imports required for tests
include("preamble_test.jl")

# main MPS and MPO modules
include("test_mps.jl")
include("test_mpo.jl")

# signals/Signals.jl
include("test_signals.jl")

# signals/SignalConverters.jl
include("test_signal_converters.jl")

# linalg/apply.jl
include("test_apply.jl")

# circuits/qft_gates.jl
include("test_qft_gates.jl")

# circuits/dt_gates.jl
include("test_dt_gates.jl")

# circuits/zt_gates.jl
include("test_zt_gates.jl")
