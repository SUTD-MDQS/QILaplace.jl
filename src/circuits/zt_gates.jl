# src/circuits/zt_gates.jl
# This module provides the ITensor gates required to implement the Zero-Trace (ZT) circuit. It imports the elementary QFT gates from qft_gates.jl and constructs them on PairedSiteMPOs that will be used in the ZT circuit.

module ZTGates
using ITensors, Printf
using ..Mpo: PairedSiteMPO
using ..QFTGates: H, P, I, Π

################################# CONTROLLED PHASE GATE MPO in 2n #####################################

# Control phase gate for qubit depth 'k' on PairedSiteMPO. In a k-qubit control gate, the QFT acts only on the copy branch while the main branch remains identity. The applicaton is similar to control_Hphase_mpo in QFTGates but only on the copy sites.
function control_Hphase_ztmps_mpo(k::Int, sites::Vector{IType}) where {IType<:Index}
    k ≥ 1 || throw(
        ArgumentError(
            "control_Hphase_ztmps_mpo: Number of qubits 'k' must be at least 1. Found k=$k",
        ),
    )
    iseven(length(sites)) || throw(
        ArgumentError(
            "control_Hphase_ztmps_mpo: Number of sites must be even. Found length(sites)=$(length(sites))",
        ),
    )
    length(sites) == 2k || throw(
        ArgumentError(
            "control_Hphase_ztmps_mpo: Number of sites must be equal to 2k. Found length(sites)=$(length(sites)), k=$k",
        ),
    )

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
    bonds_main = [Index(2; tags=@sprintf("zqft-bond-main-%d", i)) for i in 1:(k - 1)]
    bonds_copy = [Index(2; tags=@sprintf("zqft-bond-copy-%d", i)) for i in 1:k]

    # 1. Main 1 (Site 1)
    # Connects to bonds_copy[1]. Sums the branches: I * ( <1| + <2| )
    data[1] =
        I(sites_main[1]) * onehot(bonds_copy[1], 1) +
        I(sites_main[1]) * onehot(bonds_copy[1], 2)

    # 2. Copy 1 (Site 2)
    # Connects bonds_copy[1] (Left) and bonds_main[1] (Right)
    # Phase if bond=2
    # Note: We use negative angle to implement exp(-i * 2π * ...) for Z-transform/DFT
    θ = -2π / 2.0^(k)
    data[2] =
        I(sites_copy[1]) * onehot(bonds_copy[1], 1, bonds_main[1], 1) +
        P(θ, sites_copy[1]) * onehot(bonds_copy[1], 2, bonds_main[1], 2)

    # 3. Intermediate sites
    for j in 2:(k - 1)
        # Main j (Site 2j-1)
        # Connects bonds_main[j-1] (Left) and bonds_copy[j] (Right)
        # Pass through
        data[2j - 1] =
            I(sites_main[j]) * onehot(bonds_main[j - 1], 1, bonds_copy[j], 1) +
            I(sites_main[j]) * onehot(bonds_main[j - 1], 2, bonds_copy[j], 2)

        # Copy j (Site 2j)
        # Connects bonds_copy[j] (Left) and bonds_main[j] (Right)
        # Phase if bond=2
        θ = -2π / 2.0^(k-j+1)
        data[2j] =
            I(sites_copy[j]) * onehot(bonds_copy[j], 1, bonds_main[j], 1) +
            P(θ, sites_copy[j]) * onehot(bonds_copy[j], 2, bonds_main[j], 2)
    end

    # 4. Main k (Site 2k-1)
    # Connects bonds_main[k-1] (Left) and bonds_copy[k] (Right)
    # Pass through
    data[2k - 1] =
        I(sites_main[k]) * onehot(bonds_main[k - 1], 1, bonds_copy[k], 1) +
        I(sites_main[k]) * onehot(bonds_main[k - 1], 2, bonds_copy[k], 2)

    # 5. Copy k (Site 2k)
    # Connects bonds_copy[k] (Left)
    # Source of control: H then Project
    # If result 0 -> bond 1
    # If result 1 -> bond 2
    site_tmp = sim(sites_copy[k])
    H_op = replaceind!(H(sites_copy[k]), sites_copy[k]', site_tmp) # s -> tmp
    Π0_op = replaceind!(Π(0, sites_copy[k]), sites_copy[k], site_tmp) # tmp -> s'
    Π1_op = replaceind!(Π(1, sites_copy[k]), sites_copy[k], site_tmp) # tmp -> s'

    data[2k] =
        (Π0_op * H_op) * onehot(bonds_copy[k], 1) +
        (Π1_op * H_op) * onehot(bonds_copy[k], 2)

    return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
end

end # module ZTGates
