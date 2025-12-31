# src/transforms/qft_transform.jl
# This module computes the Quantum Fourier Transform (QFT) MPO that acts on a given MPS to transform it into the frequency domain. It uses the control_Hphase_mpo from QFTGates to build the full QFT MPO.

module QFTTransform
using ITensors, Printf
using ..Mpo: SingleSiteMPO
using ..Mps: SignalMPS
using ..QFTGates: control_Hphase_mpo

export build_qft_mpo

# Function to perform the zip-up algorithm that combines two control qubit MPOs and pushes the OC upwards
function zip_up_mpos(mpo1::SingleSiteMPO, mpo2::SingleSiteMPO, oc::Int)
    L1 = length(mpo1.sites)
    L2 = length(mpo2.sites)
    L1 > L2 || throw(ArgumentError("zip_up_mpos: mpo1 must be longer than mpo2. Found length(mpo1)=$L1, length(mpo2)=$L2"))
    oc == L1 || throw(ArgumentError("zip_up_mpos: Orthogonality center 'oc' must be at the end of mpo1 to do zip-up algorithm. Found oc=$oc, length(mpo1)=$L1"))

    new_data = copy(mpo1.data)
    new_bonds = copy(mpo1.bonds)

    # Initial remainder tensor (identity)
    T = ITensor(1.0)

    # Loop using i_rev from 0 to length(mpo2) - 1
    for i_rev in 0:(L2 - 1)
        idx1 = L1 - i_rev
        idx2 = L2 - i_rev

        core1 = mpo1.data[idx1]
        core2 = mpo2.data[idx2]

        # Contract cores with the remainder from the previous step to do factorisation
        core_to_factorize = core1 * core2 * T

        left_inds = Index[]
        if idx1 > 1
            push!(left_inds, mpo1.bonds[idx1-1])
        end
        if idx2 > 1
            push!(left_inds, mpo2.bonds[idx2-1])
        end

        T, V = factorize(core_to_factorize, left_inds; ortho="right", tags=@sprintf("bond-%d", idx1 - 1))

        new_data[idx1] = V
        if idx1 > 1
            new_bonds[idx1-1] = commonind(T, V)
        end
    end

    # The last remainder T becomes the first core
    new_data[L1 - L2] *= T
    
    return SingleSiteMPO(new_data, mpo1.sites, new_bonds), L1 - L2
end

# Function to perform the zip-down algorithm (SVD compression with truncation) to push the OC downwards
function zip_down_mpos(mpo::SingleSiteMPO, oc::Int; cutoff=1e-14, maxdim=1000)
    L = length(mpo.sites)
    new_data = copy(mpo.data)
    new_bonds = copy(mpo.bonds)

    # We sweep from site oc down to L-1.
    for k in oc:(L-1)
        T = new_data[k]
        
        bond_to_next = new_bonds[k]
        left_inds = uniqueinds(T, [bond_to_next])
        
        # SVD
        U, S, V = svd(T, left_inds; cutoff=cutoff, maxdim=maxdim, lefttags="bond-$(k)", righttags="bond-$(k+1)")
        new_bond = commonind(U, S)
        
        new_data[k] = U
        new_bonds[k] = new_bond
        
        # Pass remainder to next site
        remainder = S * V
        new_data[k+1] *= remainder
    end
    
    return SingleSiteMPO(new_data, mpo.sites, new_bonds), L
end

# Build the full QFT MPO for 'n' qubits using zip-up and zip-down algorithms
function build_qft_mpo(n::Int, sites::Vector{IType}; cutoff=1e-14, maxdim=1000) where {IType<:Index}
    n ≥ 1 || throw(ArgumentError("build_qft_mpo: Number of qubits 'n' must be at least 1. Found n=$n"))
    length(sites) == n || throw(ArgumentError("build_qft_mpo: Number of sites must be equal to n. Found length(sites)=$(length(sites)), n=$n"))

    if n == 1
        return control_Hphase_mpo(1, sites)
    end

    qft_mpo = control_Hphase_mpo(n, sites)
    oc = n

    # Apply the consevutive controlled phase gates
    for iter in 1:(n - 1)
        # mpo2 acts on the subset of sites (iter+1):n
        subset_sites = sites[(iter + 1):end]
        mpo2 = control_Hphase_mpo(n - iter, subset_sites)
        
        # Prepare the indices to contract inside the zip-up algo
        for (idx2, s) in enumerate(subset_sites)
            idx1 = iter + idx2
            # The main index of qft_mpo contracts with the primed index of mpo2
            s_cont = sim(s)
            qft_mpo.data[idx1] = replaceinds(qft_mpo.data[idx1], s => s_cont)
            mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
        end
        
        mpo_zipped_up, oc = zip_up_mpos(qft_mpo, mpo2, oc)
        qft_mpo, oc = zip_down_mpos(mpo_zipped_up, oc; cutoff=cutoff, maxdim=maxdim)
    end

    return qft_mpo
end

function build_qft_mpo(ψ::SignalMPS; cutoff=1e-14, maxdim=1000)
    n = length(ψ.sites)
    return build_qft_mpo(n, ψ.sites; cutoff=cutoff, maxdim=maxdim)
end

end # module QFTTransform
