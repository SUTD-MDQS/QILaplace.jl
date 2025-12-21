# src/circuits/zt_gates.jl
# This module provides the ITensor gates required to implement the Zero-Trace (ZT) circuit. It imports the elementary QFT gates from qft_gates.jl and constructs them on PairedSiteMPOs that will be used in the ZT circuit.

module ZTGates
using ITensors, Printf
using ..Mpo: PairedSiteMPO
using ..QFTGates: H, P, I, Π

################################# CONTROLLED PHASE GATE MPO in 2n #####################################

# Control phase gate for qubit depth 'k' on PairedSiteMPO. In a k-qubit control gate, the QFT acts only on the copy branch while the main branch remains identity. The applicaton is similar to control_Hphase_mpo in QFTGates but only on the copy sites.
function control_Hphase_ztmps_mpo(k::Int, sites::Vector{IType}) where {IType<:Index}
    k ≥ 1 || throw(ArgumentError("control_Hphase_ztmps_mpo: Number of qubits 'k' must be at least 1. Found k=$k"))
    iseven(length(sites)) || throw(ArgumentError("control_Hphase_ztmps_mpo: Number of sites must be even. Found length(sites)=$(length(sites))"))
    length(sites) == 2k || throw(ArgumentError("control_Hphase_ztmps_mpo: Number of sites must be equal to 2k. Found length(sites)=$(length(sites)), k=$k"))

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

    if k == 1
        # Simple Hadamard gate for single qubit on copy site and identity on main site
        data = Vector{ITensor}(undef, 2)
        sites_main = [sites[1]]
        sites_copy = [sites[2]]
        bonds_main = eltype(sites)[]
        bonds_copy = [Index(1; tags="zqft-bond-copy-1")]
        data[1] = (I(sites_main[1]) * onehot(bonds_copy[1], 1))
        data[2] = (H(sites_copy[1]) * onehot(bonds_copy[1], 1))
        return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
    end

    data = Vector{ITensor}(undef, 2k)
    sites_main = sites[1:2:end]
    sites_copy = sites[2:2:end]
    bonds_main = [Index(2, tags=@sprintf("zqft-bond-main-%d", i)) for i in 1:(k-1)]
    bonds_copy = [Index(2, tags=@sprintf("zqft-bond-copy-%d", i)) for i in 1:k]

    # Block 1 (copy[1] is Control)
    # main[1]: Identity
    data[1] = I(sites_main[1]) * onehot(bonds_copy[1], 1)
    
    # copy[1]: H + Projector.
    # Input: bonds_copy[1] (dummy 1).
    # Output: bonds_main[1] (Control State).
    site_tmp = sim(sites_copy[1])
    Hcontrol = replaceind!(H(sites_copy[1]), sites_copy[1], site_tmp)
    Π0 = replaceind!(Π(0, sites_copy[1]), sites_copy[1]', site_tmp)
    Π1 = replaceind!(Π(1, sites_copy[1]), sites_copy[1]', site_tmp)
    
    data[2] = (Hcontrol * Π0 * onehot(bonds_copy[1], 1, bonds_main[1], 1) +
               Hcontrol * Π1 * onehot(bonds_copy[1], 1, bonds_main[1], 2))

    # Intermediate sites: main = identity (pass-through), copy = controlled phase gate (reads main bond)
    for l in 2:(k-1)
        θ = 2π / 2.0^l
        # main[l]: Identity, pass signal
        data[2l-1] = (I(sites_main[l]) * onehot(bonds_main[l-1], 1, bonds_copy[l], 1) +
                      I(sites_main[l]) * onehot(bonds_main[l-1], 2, bonds_copy[l], 2))
        # copy[l]: Phase, pass signal
        data[2l] = (I(sites_copy[l]) * onehot(bonds_copy[l], 1, bonds_main[l], 1) +
                    P(θ, sites_copy[l]) * onehot(bonds_copy[l], 2, bonds_main[l], 2))
    end

    # Last site: control gate on copy site and identity on main site
    l = k
    θ = 2π / 2.0^k
    
    # main[k]: Identity, pass signal
    data[2k-1] = (I(sites_main[k]) * onehot(bonds_main[k-1], 1, bonds_copy[k], 1) + 
                  I(sites_main[k]) * onehot(bonds_main[k-1], 2, bonds_copy[k], 2))
    
    # copy[k]: Phase, end
    data[2k] = (I(sites_copy[k]) * onehot(bonds_copy[k], 1) +
                P(θ, sites_copy[k]) * onehot(bonds_copy[k], 2))

    return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
end

end # module ZTGates
