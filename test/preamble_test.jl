using Test, QILaplace
using ITensors, LinearAlgebra, Printf


# ------------------------ helper functions ------------------------

function contract_chain(data::Vector{ITensor})
    T = ITensor(1)
    for A in data
        T *= A
    end
    return T
end

function dense_from_signal_mps(ψ::QILaplace.Mps.SignalMPS)
    sites = ψ.sites
    T = contract_chain(ψ.data)
    return Array(T, sites...)
end

