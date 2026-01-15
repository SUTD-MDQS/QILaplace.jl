# src/circuits/dt_gates.jl
# This module provides the ITensor gates required to implement the Damping Transform (DT) circuit. It also constructs the controlled damping gate as an MPO that will be used in the DT circuit

module DTGates
using ITensors, Printf
using ..Mpo: PairedSiteMPO
using ..QFTGates: I, Π

export control_damping_mpo, control_damping_copy_mpo

################################ ELEMENTARY DT GATES #####################################
# Damped Hadamard gate on a qubit
dampedH(ωr::Real, site_index::IType) where {IType<:Index} = begin
    dampedHmat = 1/√2 * [1  1;
                  1 exp(-ωr / 2)]
    ITensor(dampedHmat, site_index', site_index)
end

R(ωr::Real, site_index::IType) where {IType<:Index} = begin
    Rmat = [1 0;
            0 exp(-ωr)]
    ITensor(Rmat, site_index', site_index)
end

################################# CONTROLLED DAMPING GATE MPO #####################################

# Control damping gate for qubit depth 'k'. In a k-qubit control gate, the first qubit is a control and in the rest, damped Hadamard gates apply with damping factor R: ωr/2^2 to 2^k in the main branch and I in the copy branch
function control_damping_mpo(n::Int, k::Int, ωr::Real, sites::Vector{IType}) where {IType<:Index}
    k ≥ 1 || throw(ArgumentError("control_damping_mpo: Number of qubits 'k' must be at least 1. Found k=$k"))
    iseven(length(sites)) || throw(ArgumentError("control_damping_mpo: Number of sites must be even. Found length(sites)=$(length(sites))"))
    length(sites) == 2k || throw(ArgumentError("control_damping_mpo: Number of sites must be equal to 2k. Found length(sites)=$(length(sites)), k=$k"))

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
        # Simple damped Hadamard gate for main qubit and identity on copy qubit
        data = Vector{ITensor}(undef, 2)
        sites_main = [sites[1]]
        sites_copy = [sites[2]]
        bonds_main = eltype(sites)[]
        bonds_copy = [Index(1, tags=@sprintf("dt-bond-copy-%d", i)) for i in 1:k]
        data[1] = dampedH(ωr, sites[1]) * onehot(bonds_copy[1], 1)
        data[2] = I(sites[2]) * onehot(bonds_copy[1], 1)
        return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
    end

    data = Vector{ITensor}(undef, 2k)
    sites_main = sites[1:2:end]
    sites_copy = sites[2:2:end]
    bonds_main = [Index(2, tags=@sprintf("dt-bond-main-%d", i)) for i in 1:(k-1)]
    bonds_copy = [Index(2, tags=@sprintf("dt-bond-copy-%d", i)) for i in 1:k]
    
    # Intermediate sites: controlled R and identity gates
    for l in 1:(k-1)
        R_factor = ωr * 2.0^(l - k - 1)
        # main site tensor: controlled R gate
        data[2l-1] = begin
            if l == 1
                # Onehot-encode just for the copy bond
                (I(sites_main[l]) * onehot(bonds_copy[l], 1) +
                 R(R_factor, sites_main[l]) * onehot(bonds_copy[l], 2))
            else
                # Onehot-encode for both main and copy bonds
                (I(sites_main[l]) * onehot(bonds_main[l-1], 1, bonds_copy[l], 1) +
                      R(R_factor, sites_main[l]) * onehot(bonds_main[l - 1], 2, bonds_copy[l], 2))
            end
        end
3
        # copy site tensor: identity gate
        data[2l] = (I(sites_copy[l]) * onehot(bonds_copy[l], 1, bonds_main[l], 1) + 
                    I(sites_copy[l]) * onehot(bonds_copy[l], 2, bonds_main[l], 2))
    end

    # Last site: copy gate and closing bonds
    # Control Logic: |0> -> Bond 1, |1> -> Bond 2
    Π0 = Π(0, sites_main[k])
    Π1 = Π(1, sites_main[k])
    Hdamped = dampedH(ωr, sites_main[k])
    site_tmp = sim(sites_main[k])
    replaceind!(Π0, sites_main[k], site_tmp)
    replaceind!(Π1, sites_main[k], site_tmp)
    replaceind!(Hdamped, sites_main[k]', site_tmp)

    data[2k-1] = ((Π0 * Hdamped * onehot(bonds_main[k-1], 1, bonds_copy[k], 1) + 
                    Π1 * Hdamped * onehot(bonds_main[k-1], 2, bonds_copy[k], 2)))

    # Last site: copy gate
    # This is the end of the block, so no right bonds (bonds_main[k] does not exist)
    data[2k] = (I(sites_copy[k])  * onehot(bonds_copy[k], 1) +
                I(sites_copy[k])  * onehot(bonds_copy[k], 2))

    return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
end


function control_damping_copy_mpo(n::Int, k::Int, ωr::Real, sites::Vector{IType}) where {IType<:Index}
    k ≥ 1 || throw(ArgumentError("control_damping_copy_mpo: Number of qubits 'k' must be at least 1. Found k=$k"))
    iseven(length(sites)) || throw(ArgumentError("control_damping_copy_mpo: Number of sites must be even. Found length(sites)=$(length(sites))"))
    
    # The MPO acts on qubits k through n.
    # Number of qubits involved: L = n - k + 1
    L = n - k + 1
    length(sites) == 2L || throw(ArgumentError("control_damping_copy_mpo: Number of sites must be equal to 2*(n-k+1). Found length(sites)=$(length(sites)), n=$n, k=$k"))

    # one-hot tensor embedding for bond indices
    onehot(b::Index, j::Int) = begin
        T = ITensor(b)
        T[b => j] = 1
        return T
    end
    onehot(bL::Index, r::Int, bR::Index, c::Int) = begin
        T = ITensor(bL, bR)
        T[bL => r, bR => c] = 1
        return T
    end

    data = Vector{ITensor}(undef, 2L)
    sites_main = sites[1:2:end]
    sites_copy = sites[2:2:end]
    
    if L == 1
        # k=n: Just identity on main[1] and copy[1] (relative indices)
        # This corresponds to global qubit n.
        bonds_copy = [Index(1, tags=@sprintf("dt-bond-copy-%d", 1))]
        data[1] = I(sites_main[1]) * onehot(bonds_copy[1], 1)
        data[2] = I(sites_copy[1]) * onehot(bonds_copy[1], 1)
        return PairedSiteMPO(data, sites_main, sites_copy, eltype(sites)[], bonds_copy)
    end
    
    bonds_main = [Index(2, tags=@sprintf("dt-bond-main-%d", i)) for i in 1:(L-1)]
    bonds_copy = [Index(2, tags=@sprintf("dt-bond-copy-%d", i)) for i in 1:L]
    
    # main[1]: Identity
    data[1] = I(sites_main[1]) * onehot(bonds_copy[1], 1)
    
    # copy[1]: Control Projector
    # Control Logic: |0X0| -> Bond 1, |1X1| -> Bond 2
    data[2] = (Π(0, sites_copy[1]) * onehot(bonds_copy[1], 1, bonds_main[1], 1) +
               Π(1, sites_copy[1]) * onehot(bonds_copy[1], 1, bonds_main[1], 2))
    
    # Blocks 2 to L-1 (Global k+1 to n-1)
    for j in 2:(L-1)
        # R-factor for target j (Global k+j-1)
        # θ = ωr * 2^(j-2)
        R_factor = ωr * 2.0^(j - 2)
        
        # main[j]: Controlled R
        data[2j-1] = (I(sites_main[j]) * onehot(bonds_main[j-1], 1, bonds_copy[j], 1) +
                      R(R_factor, sites_main[j]) * onehot(bonds_main[j-1], 2, bonds_copy[j], 2))
        # copy[j]: Identity
        data[2j] = (I(sites_copy[j]) * onehot(bonds_copy[j], 1, bonds_main[j], 1) +
                    I(sites_copy[j]) * onehot(bonds_copy[j], 2, bonds_main[j], 2))
    end
    
    # Last Block L (Global n)
    j = L
    R_factor = ωr * 2.0^(j - 2)
    
    # main[L]: Controlled R
    data[2L-1] = (I(sites_main[L]) * onehot(bonds_main[L-1], 1, bonds_copy[L], 1) +
                  R(R_factor, sites_main[L]) * onehot(bonds_main[L-1], 2, bonds_copy[L], 2))
    
    # copy[L]: Identity
    data[2L] = (I(sites_copy[L]) * onehot(bonds_copy[L], 1) +
                I(sites_copy[L]) * onehot(bonds_copy[L], 2))
    
    return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
end

end # module DTGates
