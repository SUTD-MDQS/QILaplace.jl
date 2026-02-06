# src/QILaplace.jl
module QILaplace
__precompile__(true)

using ITensors, Random, Printf

# Mps.jl
include("mps.jl")
using .Mps:
    SignalMPS,
    zTMPS,
    nsite,
    siteindices,
    bondindices,
    canonicalize!,
    compress!,
    update_site!,
    update_bond!,
    coefficient
export SignalMPS,
    zTMPS,
    nsite,
    siteindices,
    bondindices,
    canonicalize!,
    compress!,
    update_site!,
    update_bond!,
    coefficient

# Mpo.jl
include("mpo.jl")
using .Mpo: SingleSiteMPO, PairedSiteMPO
export SingleSiteMPO, PairedSiteMPO

# linalg/apply.jl
include("linalg/apply.jl")
using .ApplyMPO: apply

# signals/Signals.jl
include("signals/Signals.jl")
using .Signals: generate_signal
export generate_signal

# linalg/rsvd.jl
include("linalg/rsvd.jl")
using .RSVD: rsvd

# signals/SignalConverters.jl
include("signals/SignalConverters.jl")
using .SignalConverters: signal_mps, signal_ztmps
export signal_mps, signal_ztmps

# circuits/qft_gates.jl
include("circuits/qft_gates.jl")
using .QFTGates

# circuits/dt_gates.jl
include("circuits/dt_gates.jl")
using .DTGates

# circuits/zt_gates.jl
include("circuits/zt_gates.jl")
using .ZTGates

# transforms/qft_transformer.jl
include("transforms/qft_transformer.jl")
using .QFTTransform: build_qft_mpo
export build_qft_mpo

# transforms/dt_transformer.jl
include("transforms/dt_transformer.jl")
using .DTTransform: build_dt_mpo
export build_dt_mpo

# transforms/zt_transformer.jl
include("transforms/zt_transformer.jl")
using .ZTTransformer: build_zt_mpo
export build_zt_mpo

# Algorithms.jl

# Precompile
__init__() = nothing

end # module QILaplace
