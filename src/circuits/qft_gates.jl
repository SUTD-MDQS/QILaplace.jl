# src/ciruits/qft_gates.jl
# This module provides the ITensor gates required to implement the Quantum Fourier Transform (QFT) circuit. It also constructs the controlled phase rotation gate as an MPO that will be used in the QFT circuit.

module QFTGates

using ITensors, Printf
using ..Mpo: SingleSiteMPO

export control_Hphase_mpo

################################ ELEMENTARY QFT GATES #####################################

# The identity operation on a qubit
I(site_index::IType) where {IType<:Index} = delta(site_index', site_index)

# The Hadamard gate on a qubit
H(site_index::IType) where {IType<:Index} = begin
    Hmat = (1/sqrt(2)) * [1  1;
                        1 -1]
    ITensor(Hmat, site_index', site_index)
end

# The phase rotation gate between two qubits
P(θ, site_index::IType) where {IType<:Index} = begin
    Pmat = [1 0;
            0 exp(im * θ)]
    ITensor(Pmat, site_index', site_index)
end

Π(i::Int, site_index::IType) where {IType<:Index} = begin
    i < dim(site_index) || throw(ArgumentError("Π: Index dimension is less than $i"))
    Projector = ITensor(site_index', site_index)
    i_idx = i + 1 # ITensors are 1-indexed
    Projector[i_idx, i_idx] = 1.0
    return Projector
end

################################# CONTROLLED PHASE GATE MPO #####################################

# Control gate for qubit depth 'k'. In a k-qubit control gate, the first qubit is a control and in the rest, phase gates apply rotation starting from 2π/2^2 to 2π/2^k
function control_Hphase_mpo(k::Int, sites::Vector{IType}) where {IType<:Index}
    k ≥ 1 || throw(ArgumentError("control_Hphase_mpo: Number of qubits 'k' must be at least 1. Found k=$k"))
    length(sites) == k || throw(ArgumentError("control_Hphase_mpo: Number of sites must be equal to k. Found length(sites)=$(length(sites)), k=$k"))
    if k == 1
        # Simple Hadamard gate for single qubit
        bonds = eltype(sites)[]
        return SingleSiteMPO([H(sites[1])], sites, bonds)
    end

    # one-hot tensor embedding for bond indices
    onehot(b::Index, k::Int) = begin
        T = ITensor(b)
        T[b => k] = 1
        return T
    end
    onehot(bL::Index, r::Int, bR::Index, c::Int) = begin
        T = ITensor(bL, bR)
        T[bL => r, bR => c] = 1
        return T
    end

    data = Vector{ITensor}(undef, k)
    bonds =  [Index(2; tags=@sprintf("qft-bond-%d", i)) for i in 1:(k - 1)]

    # First site: control projector and bond to next site
    site_tmp = sim(sites[1]) # To ensure unique tags for contraction
    Hcontrol = replaceind!(H(sites[1]), sites[1], site_tmp)
    Π0 = replaceind!(Π(0, sites[1]), sites[1]', site_tmp)
    Π1 = replaceind!(Π(1, sites[1]), sites[1]', site_tmp)
    data[1] = (Hcontrol * Π0 * onehot(bonds[1], 1) +
              Hcontrol * Π1 * onehot(bonds[1], 2))
    
    # Intermediate sites: controlled phase rotations
    for l in 2:(k - 1)
        θ = 2π / 2.0^(l)
        data[l] = (I(sites[l]) * onehot(bonds[l - 1], 1, bonds[l], 1) +
                  P(θ, sites[l]) * onehot(bonds[l - 1], 2, bonds[l], 2))
    end
    # Last site: phase rotation and closing bond
    θ = 2π / 2.0^(k)
    data[k] = (I(sites[k]) * onehot(bonds[k - 1], 1) +
              P(θ, sites[k]) * onehot(bonds[k - 1], 2))

    return SingleSiteMPO(data, sites, bonds)
end

end # module QFTGates
