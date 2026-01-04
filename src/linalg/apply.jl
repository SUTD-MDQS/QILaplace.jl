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
        # Define the site to contract
        site_tmp = sim(ψ.sites[i])
        # Prepare MPO tensor by replacing site indices
        T_mpo = copy(W.data[i])
        T_mps = copy(ψ.data[i])
        replaceind!(T_mpo, W.sites[i]', site_tmp)
        replaceind!(T_mps, ψ.sites[i], site_tmp)
        # Contract MPO tensor with MPS tensor
        new_data[i] = T_mpo * T_mps
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

"""
Apply MPO to MPO (contract matching sites, combine bonds).
Equivalent to mpo1 * mpo2 for PairedSiteMPOs. This function:
 - aligns MPOs based on site indices (using `id` to match)
 - contracts tensors where sites match
 - combines bonds where both MPOs have bonds
 - propagates single bonds where only one MPO is present
"""
function apply(mpo1::PairedSiteMPO, mpo2::PairedSiteMPO)
    # Map site ID to index in MPO
    map1 = Dict(id(s) => i for (i, s) in enumerate(mpo1.sites_main))
    map2 = Dict(id(s) => i for (i, s) in enumerate(mpo2.sites_main))
    
    # Determine union of sites, preserving order from mpo1 then mpo2
    all_sites = Vector{eltype(mpo1.sites_main)}()
    seen_ids = Set{ITensors.IDType}()
    
    for s in mpo1.sites_main
        push!(all_sites, s)
        push!(seen_ids, id(s))
    end
    
    for s in mpo2.sites_main
        if !(id(s) in seen_ids)
            push!(all_sites, s)
            push!(seen_ids, id(s))
        end
    end
    
    N = length(all_sites)
    new_data = Vector{ITensor}(undef, 2N)
    new_sites_main = all_sites
    new_sites_copy = Vector{eltype(mpo1.sites_copy)}(undef, N)
    
    # Populate new_sites_copy
    for (i, s) in enumerate(all_sites)
        if haskey(map1, id(s))
            new_sites_copy[i] = mpo1.sites_copy[map1[id(s)]]
        elseif haskey(map2, id(s))
            new_sites_copy[i] = mpo2.sites_copy[map2[id(s)]]
        else
            error("Site not found in either MPO")
        end
    end
    
    new_bonds_main = Vector{eltype(mpo1.sites_main)}(undef, N-1)
    new_bonds_copy = Vector{eltype(mpo1.sites_main)}(undef, N)
    
    C_main_prev = nothing # Combiner from previous block (main bond)
    
    for i in 1:N
        s = all_sites[i]
        idx1 = get(map1, id(s), nothing)
        idx2 = get(map2, id(s), nothing)
        
        # --- MAIN TENSOR ---
        M1 = (idx1 !== nothing) ? copy(mpo1.data[2*idx1-1]) : nothing
        M2 = (idx2 !== nothing) ? copy(mpo2.data[2*idx2-1]) : nothing
        
        if M1 !== nothing && M2 !== nothing
            s_tmp = sim(s)
            replaceind!(M1, s, s_tmp)
            replaceind!(M2, s', s_tmp)
            M_new = M1 * M2
        elseif M1 !== nothing
            M_new = M1
        elseif M2 !== nothing
            M_new = M2
        else
            error("Site not found in either MPO")
        end
        
        if C_main_prev !== nothing
            M_new = M_new * C_main_prev
        end
        
        # Combine internal copy bonds
        b_c1 = (idx1 !== nothing) ? mpo1.bonds_copy[idx1] : nothing
        b_c2 = (idx2 !== nothing) ? mpo2.bonds_copy[idx2] : nothing
        
        C_copy = nothing
        if b_c1 !== nothing && b_c2 !== nothing
            C_copy = combiner(b_c1, b_c2; tags = "bond-copy-$(i)")
            M_new = M_new * C_copy
            new_bonds_copy[i] = combinedind(C_copy)
        elseif b_c1 !== nothing
            new_bonds_copy[i] = b_c1
        elseif b_c2 !== nothing
            new_bonds_copy[i] = b_c2
        else
             error("Missing copy bond at site $i")
        end
        
        new_data[2i-1] = M_new
        
        # --- COPY TENSOR ---
        C1 = (idx1 !== nothing) ? copy(mpo1.data[2*idx1]) : nothing
        C2 = (idx2 !== nothing) ? copy(mpo2.data[2*idx2]) : nothing
        
        s_c = new_sites_copy[i]
        if C1 !== nothing && C2 !== nothing
            s_tmp = sim(s_c)
            replaceind!(C1, s_c, s_tmp)
            replaceind!(C2, s_c', s_tmp)
            C_new = C1 * C2
        elseif C1 !== nothing
            C_new = C1
        elseif C2 !== nothing
            C_new = C2
        else
            error("Copy site tensor missing")
        end
        
        if C_copy !== nothing
            C_new = C_new * dag(C_copy)
        end
        
        # Combine main bonds (to next block)
        if i < N
            b_m1 = (idx1 !== nothing && idx1 < length(mpo1.sites_main)) ? mpo1.bonds_main[idx1] : nothing
            b_m2 = (idx2 !== nothing && idx2 < length(mpo2.sites_main)) ? mpo2.bonds_main[idx2] : nothing
            
            if b_m1 !== nothing && b_m2 !== nothing
                C_main_prev = combiner(b_m1, b_m2; tags = "bond-main-$(i)")
                C_new = C_new * C_main_prev
                new_bonds_main[i] = combinedind(C_main_prev)
            elseif b_m1 !== nothing
                new_bonds_main[i] = b_m1
                C_main_prev = nothing
            elseif b_m2 !== nothing
                new_bonds_main[i] = b_m2
                C_main_prev = nothing
            else
                # Disconnected chain or gap?
                # If we are here, i < N, so there is a next site.
                # If neither has a bond, we can't connect.
                error("Disconnected MPO chain at site $i")
            end
        else
            C_main_prev = nothing
        end
        
        new_data[2i] = C_new
    end
    
    return PairedSiteMPO(new_data, new_sites_main, new_sites_copy, new_bonds_main, new_bonds_copy)
end

end # module ApplyMPO
