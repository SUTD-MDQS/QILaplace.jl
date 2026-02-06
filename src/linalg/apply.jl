# src/linalg/apply.jl
# This module provides the function that applies an MPO to an MPS. The output is an MPS with updated site indices.

module ApplyMPO

using ITensors, Printf
import Base: *
using ..Mps
using ..Mps: SignalMPS, zTMPS
using ..Mpo: SingleSiteMPO, PairedSiteMPO

export apply

function _as_single_site_mpo(W::PairedSiteMPO)
    n = length(W.sites_main)
    sites = Vector{eltype(W.sites_main)}(undef, 2n)
    bonds = Vector{eltype(W.bonds_main)}(undef, 2n - 1)

    for i in 1:n
        sites[2i - 1] = W.sites_main[i]
        sites[2i] = W.sites_copy[i]

        bonds[2i - 1] = W.bonds_copy[i]
        if i < n
            bonds[2i] = W.bonds_main[i]
        end
    end

    return SingleSiteMPO(W.data, sites, bonds)
end

function _paired_from_single(W::SingleSiteMPO)
    iseven(length(W)) || throw(
        ArgumentError(
            "_paired_from_single: length must be even to split into PairedSiteMPO."
        ),
    )
    n = length(W) ÷ 2

    sites_main = Vector{eltype(W.sites)}(undef, n)
    sites_copy = Vector{eltype(W.sites)}(undef, n)
    bonds_main = Vector{eltype(W.bonds)}(undef, max(n - 1, 0))
    bonds_copy = Vector{eltype(W.bonds)}(undef, n)

    for i in 1:n
        sites_main[i] = W.sites[2i - 1]
        sites_copy[i] = W.sites[2i]
        bonds_copy[i] = W.bonds[2i - 1]
        if i < n
            bonds_main[i] = W.bonds[2i]
        end
    end

    new_data = [copy(t) for t in W.data]
    return PairedSiteMPO(new_data, sites_main, sites_copy, bonds_main, bonds_copy)
end

# Function to apply an MPO to an MPS
function apply(W::SingleSiteMPO, ψ::SignalMPS; kwargs...)
    length(W) == length(ψ) || throw(
        ArgumentError(
            "apply: MPO and MPS must have the same number of sites. Found length(W)=$(length(W)), length(psi)=$(length(ψ))",
        ),
    )
    W.sites == ψ.sites || throw(
        ArgumentError(
            "apply: MPO and MPS must have the same site indices. Found W.sites=$(W.sites), psi.sites=$(ψ.sites)",
        ),
    )

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
    for i in 1:(n - 1)
        # Create combiner for the MPO bond and MPS bond
        combining_inds = (W.bonds[i], ψ.bonds[i])
        combiner_tensor = combiner(combining_inds...)
        combined_ind = combinedind(combiner_tensor)
        new_bond = sim(combined_ind; tags=@sprintf("bond-%d", i))
        replaceind!(combiner_tensor, combined_ind, new_bond)

        # Apply combiner to the right side of site i
        new_data[i] = new_data[i] * combiner_tensor
        new_bonds[i] = new_bond

        # Apply dag combiner to the left side of site i+1
        new_data[i + 1] = dag(combiner_tensor) * new_data[i + 1]
    end

    return SignalMPS(new_data, ψ.sites, new_bonds)
end

function apply(W1::SingleSiteMPO, W2::SingleSiteMPO; kwargs...)
    n1 = length(W1)
    n2 = length(W2)

    # 1. Inspect window
    start1 = findfirst(s -> s in W2.sites, W1.sites)
    start1 === nothing && throw(ArgumentError("apply: No matching sites found"))
    start2 = findfirst(s -> s == W1.sites[start1], W2.sites)

    match_len = 0
    while start1 + match_len <= n1 &&
              start2 + match_len <= n2 &&
              W1.sites[start1 + match_len] == W2.sites[start2 + match_len]
        match_len += 1
    end

    # 2. Determine base MPO (longer or W1 if equal)
    if n1 >= n2
        W_base, W_other = W1, W2
        base_start, other_start = start1, start2
    else
        W_base, W_other = W2, W1
        base_start, other_start = start2, start1
    end

    # 3. Initialize result from base MPO (preserves non-overlapping regions)
    n_out = length(W_base)
    new_data = copy(W_base.data)
    new_sites = copy(W_base.sites)
    new_bonds = copy(W_base.bonds)

    # 4. Process matching window
    prev_comb = nothing

    for i in 0:(match_len - 1)
        idx1 = start1 + i
        idx2 = start2 + i
        base_idx = base_start + i

        # Contract tensors: W1 site output (unprimed) -> W2 site input (primed)
        site_tmp = sim(W1.sites[idx1])
        T1 = copy(W1.data[idx1])
        T2 = copy(W2.data[idx2])

        replaceind!(T1, W1.sites[idx1], site_tmp)
        replaceind!(T2, W2.sites[idx2]', site_tmp)

        T = T1 * T2

        # Apply previous combiner (left side of T)
        if prev_comb !== nothing
            T = T * dag(prev_comb)
        end

        # Create/Apply current combiner (right side of T) - if not last in window
        if i < match_len - 1
            b1 = W1.bonds[idx1]
            b2 = W2.bonds[idx2]

            comb = combiner(b1, b2)
            cind = combinedind(comb)
            new_bond = sim(cind; tags=@sprintf("bond-%d", base_idx))
            replaceind!(comb, cind, new_bond)

            T = T * comb
            new_bonds[base_idx] = new_bond
            prev_comb = comb
        else
            prev_comb = nothing
        end

        new_data[base_idx] = T
    end

    return SingleSiteMPO(new_data, new_sites, new_bonds)
end

function apply(W::PairedSiteMPO, ψ::zTMPS; kwargs...)
    length(W.data) == 2 * length(ψ.sites_main) ||
        throw(ArgumentError("apply: MPO and MPS must have compatible sizes."))

    # Convert zTMPS to 2n SignalMPS
    ψ_2n = Mps._as_signal_2n(ψ)

    # Convert PairedSiteMPO to SingleSiteMPO (2n sites)
    W_single = _as_single_site_mpo(W)

    # Apply using the SingleSiteMPO implementation
    ψ_out_2n = apply(W_single, ψ_2n; kwargs...)

    # Convert back to zTMPS
    return Mps._writeback_signal_2n(ψ_out_2n)
end

function apply(mpo1::PairedSiteMPO, mpo2::PairedSiteMPO; kwargs...)
    # Convert to SingleSiteMPO representation (2n-site form)
    mpo1_single = _as_single_site_mpo(mpo1)
    mpo2_single = _as_single_site_mpo(mpo2)

    # Apply using SingleSiteMPO logic (handles unequal lengths)
    combined_single = apply(mpo1_single, mpo2_single; kwargs...)

    # Convert back to PairedSiteMPO
    return _paired_from_single(combined_single)
end

# Convenience operator overloads
*(W::SingleSiteMPO, ψ::SignalMPS) = apply(W, ψ)
*(W1::SingleSiteMPO, W2::SingleSiteMPO) = apply(W1, W2)
*(W::PairedSiteMPO, ψ::zTMPS) = apply(W, ψ)
*(W1::PairedSiteMPO, W2::PairedSiteMPO) = apply(W1, W2)

end # module ApplyMPO
