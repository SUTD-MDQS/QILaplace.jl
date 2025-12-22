# src/QILaplace.jl
module QILaplace
__precompile__(true)

using ITensors, Random, Printf

# Mps.jl
include("Mps.jl")
using .Mps: SignalMPS, zTMPS,
            nsite, siteindices, bondindices,
            canonicalize!, compress!,
            update_site!, update_bond!
export SignalMPS, zTMPS,
        nsite, siteindices, bondindices,
        canonicalize!, compress!,
        update_site!, update_bond!

# Mpo.jl
include("Mpo.jl")
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

# Transforms.jl

# Algorithms.jl

# Precompile
__init__() = nothing

end # module QILaplace
