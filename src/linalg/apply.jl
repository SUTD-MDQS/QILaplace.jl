# src/linalg/apply.jl
# This module provides the function that applies an MPO to an MPS. The output is an MPS with updated site indices.

module ApplyMPO

using ITensors, Printf
using ..Mps
using ..Mps: SignalMPS, zTMPS
using ..Mpo: SingleSiteMPO, PairedSiteMPO

export apply

# Function to apply an MPO to an MPS
function apply(W::SingleSiteMPO, ψ::SignalMPS; kwargs...)
    length(W) == length(ψ) || throw(ArgumentError("apply: MPO and MPS must have the same number of sites. Found length(W)=$(length(W)), length(psi)=$(length(ψ))"))
    W.sites == ψ.sites || throw(ArgumentError("apply: MPO and MPS must have the same site indices. Found W.sites=$(W.sites), psi.sites=$(ψ.sites)"))

    n = length(ψ)
    IndexType = eltype(ψ.sites)
    new_data = Vector{ITensor}(undef, n)
    new_bonds = Vector{IndexType}(undef, n - 1)

    for i in 1:n
        # Contract MPO tensor with MPS tensor
        new_data[i] = W.data[i] * ψ.data[i]
        # Replace primed site index with unprimed one
        replaceind!(new_data[i], ψ.sites[i]', ψ.sites[i])
    end

    # Combine the MPO and MPS bond indices at each bond
    for i in 1:n-1
        # Create combiner for the MPO bond and MPS bond
        combining_inds = (W.bonds[i], ψ.bonds[i])
        combiner_tensor = combiner(combining_inds...)
        combined_ind = combinedind(combiner_tensor)
        new_bond = sim(combined_ind, tags=@sprintf("bond-%d", i))
        replaceind!(combiner_tensor, combined_ind, new_bond)

        # Apply combiner to the right side of site i
        new_data[i] = new_data[i] * combiner_tensor
        new_bonds[i] = new_bond

        # Apply dag combiner to the left side of site i+1
        new_data[i+1] = dag(combiner_tensor) * new_data[i+1]
    end

    return SignalMPS(new_data, ψ.sites, new_bonds)
end

function _as_single_site_mpo(W::PairedSiteMPO)
    n = length(W.sites_main)
    sites = Vector{eltype(W.sites_main)}(undef, 2n)
    bonds = Vector{eltype(W.bonds_main)}(undef, 2n - 1)
    
    for i in 1:n
        sites[2i-1] = W.sites_main[i]
        sites[2i] = W.sites_copy[i]
        
        bonds[2i-1] = W.bonds_copy[i]
        if i < n
            bonds[2i] = W.bonds_main[i]
        end
    end
    
    return SingleSiteMPO(W.data, sites, bonds)
end

function apply(W::PairedSiteMPO, ψ::zTMPS; kwargs...)
    length(W.data) == 2 * length(ψ.sites_main) || throw(ArgumentError("apply: MPO and MPS must have compatible sizes."))
    
    # Convert zTMPS to 2n SignalMPS
    ψ_2n = Mps._as_signal_2n(ψ)
    
    # Convert PairedSiteMPO to SingleSiteMPO (2n sites)
    W_single = _as_single_site_mpo(W)
    
    # Apply using the SingleSiteMPO implementation
    ψ_out_2n = apply(W_single, ψ_2n; kwargs...)
    
    # Convert back to zTMPS
    return Mps._writeback_signal_2n(ψ_out_2n)
end

end # module ApplyMPO
