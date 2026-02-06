# src/transforms/dt_transform.jl
# This module builds the Damping Transform MPO that acts on a given MPS to apply the damping transform circuit. 
# It uses the dt_gates from DTGates to build the full DT MPO as a PairedSiteMPO that acts on zTMPS representations of signals.

# The DT transform has two parts:
# - Part 1: control_damping_mpo blocks with control on main site m, targets on main sites 1..m-1
# - Part 2: control_damping_copy_mpo blocks with control on copy site m, targets on main sites m+1..n

# The algorithm composes these blocks using zip-up and zip-down compression.

module DTTransform
using ITensors, Printf, LinearAlgebra
using ..Mpo: PairedSiteMPO
using ..Mps: zTMPS
using ..DTGates: control_damping_mpo, control_damping_copy_mpo

export build_dt_mpo

# Function to do zip-to-combine of two mpos
function zip_to_combine_mpos(mpo1::PairedSiteMPO, mpo2::PairedSiteMPO, oc::Int)
    L1 = length(mpo1.sites_main)
    L2 = length(mpo2.sites_main)

    @assert L1 >= L2 "zip_to_combine_mpos: mpo1 must be longer than mpo2. Found length(mpo1)=$L1, length(mpo2)=$L2"

    new_data = copy(mpo1.data)
    new_sites_main = copy(mpo1.sites_main)
    new_sites_copy = copy(mpo1.sites_copy)
    new_bonds_main = copy(mpo1.bonds_main)
    new_bonds_copy = copy(mpo1.bonds_copy)

    # Initial remainder tensor (identity)
    T = ITensor(1.0)

    # Decide which way to go in order to combine the sites
    direction::String = ""

    if mpo1.sites_main[1] == mpo2.sites_main[1] # Zip-down (aligned at start)
        direction = "down"
        # Loop over all shared tensors (main and copy pairs)
        for k in 1:(2 * L2)
            idx1 = k
            idx2 = k

            core1 = mpo1.data[idx1] # Indices: (s, s', l1, r1)
            core2 = mpo2.data[idx2] # Indices: (s, s', l2, r2)

            # Determine physical indices
            is_main = isodd(k)
            site_idx = is_main ? (k+1)÷2 : k÷2
            s = is_main ? mpo1.sites_main[site_idx] : mpo1.sites_copy[site_idx]
            s_prime = s'

            # We want to compute mpo2 * mpo1 (mpo1 acts first, then mpo2)
            # Contract mpo1's output (s') with mpo2's input (s)
            link = sim(s)
            core1_mod = replaceind(core1, s => link)
            core2_mod = replaceind(core2, s_prime => link)

            # Contract
            core = core2_mod * core1_mod * T

            # Factorize
            # Left indices: physical indices + bond from previous
            left_inds = Index[s, s_prime]
            if k > 1
                # Find the bond connecting to previous tensor
                prev_tensor = new_data[k - 1]
                common_bonds = commoninds(prev_tensor, core)
                append!(left_inds, common_bonds)
            end

            U, S, V = svd(core, left_inds...; cutoff=1e-14, maxdim=1000)

            new_data[idx1] = U
            T = S * V

            # Update bonds
            bond = commonind(U, S)
            if is_main
                # main[i] -> copy[i]
                new_bonds_copy[site_idx] = bond
            else
                # copy[i] -> main[i+1]
                if site_idx < L1
                    new_bonds_main[site_idx] = bond
                end
            end
        end

        # Absorb remainder into the next tensor of mpo1 (if it exists)
        if L1 > L2
            new_data[2 * L2 + 1] *= T
        else
            new_data[2 * L2] *= T
        end

    elseif mpo1.sites_main[end] == mpo2.sites_main[end] # Zip-up (aligned at end)
        direction = "up"
        # Loop backwards over shared tensors
        for k in 0:(2 * L2 - 1)
            idx1 = 2*L1 - k
            idx2 = 2*L2 - k

            core1 = mpo1.data[idx1]
            core2 = mpo2.data[idx2]

            # Determine physical indices
            is_main = isodd(idx1)
            site_idx = is_main ? (idx1+1)÷2 : idx1÷2
            s = is_main ? mpo1.sites_main[site_idx] : mpo1.sites_copy[site_idx]
            s_prime = s'

            # Compute mpo2 * mpo1
            link = sim(s)
            core1_mod = replaceind(core1, s => link)
            core2_mod = replaceind(core2, s_prime => link)

            core = core2_mod * core1_mod * T

            # Factorize
            # Right indices: physical indices + bond from right (next processed, which is idx+1)
            right_inds = Index[s, s_prime]
            if k > 0
                prev_tensor = new_data[idx1 + 1] # "Previous" in processing order (right neighbor)
                common_bonds = commoninds(prev_tensor, core)
                append!(right_inds, common_bonds)
            end

            # Factorize core -> T * V (T is remainder to left, V is site tensor)
            # We want V to have right_inds.
            rem_inds = uniqueinds(core, right_inds)
            U, S, V = svd(core, rem_inds...; cutoff=1e-14, maxdim=1000)

            new_data[idx1] = V
            T = U * S

            # Update bonds
            bond = commonind(S, V)
            if !is_main # copy[i]
                # main[i] <- copy[i]
                new_bonds_copy[site_idx] = bond
            else # main[i]
                # copy[i-1] <- main[i]
                if site_idx > 1
                    new_bonds_main[site_idx - 1] = bond
                end
            end
        end

        # Absorb remainder into the left neighbor of mpo1 (if it exists)
        if L1 > L2
            new_data[2 * L1 - 2 * L2] *= T
        else
            new_data[2 * L1 - 2 * L2 + 1] *= T
        end

    else
        throw(ArgumentError("zip_to_combine_mpos: Unable to determine zip direction."))
    end

    oc = L1 - L2 # This is just a placeholder, OC tracking needs to be consistent with usage
    return PairedSiteMPO(
        new_data, new_sites_main, new_sites_copy, new_bonds_main, new_bonds_copy
    ),
    oc,
    direction
end

# Function to compress the zipped MPO starting from the given oc
function zip_to_compress_mpo(
    mpo::PairedSiteMPO, oc::Int, direction::String; cutoff=1e-14, maxdim=1000
)
    L = length(mpo.data) # Total tensors (2n)
    new_data = copy(mpo.data)
    new_bonds_main = copy(mpo.bonds_main)
    new_bonds_copy = copy(mpo.bonds_copy)

    if direction == "down"
        # Move OC to the end (sweep 1 -> L)
        for i in 1:(L - 1)
            # Contract i and i+1
            core = new_data[i] * new_data[i + 1]

            # SVD
            # Left indices: unique indices of i (excluding bond to i+1)
            left_inds = uniqueinds(new_data[i], new_data[i + 1])
            U, S, V = svd(core, left_inds...; cutoff=cutoff, maxdim=maxdim)

            new_data[i] = U
            new_data[i + 1] = S * V

            # Update bond
            bond = commonind(U, S)
            # i is index in data (1..2n)
            # If i is odd (main[m]), bond is bonds_copy[m]
            # If i is even (copy[c]), bond is bonds_main[c]
            if isodd(i)
                m = (i+1)÷2
                new_bonds_copy[m] = bond
            else
                c = i÷2
                if c <= length(new_bonds_main)
                    new_bonds_main[c] = bond
                end
            end
        end
        oc = L

    elseif direction == "up"
        # Move OC to the start (sweep L -> 1)
        for i in L:-1:2
            # Contract i and i-1
            core = new_data[i] * new_data[i - 1]

            # SVD
            # Right indices: unique indices of i (excluding bond to i-1)
            right_inds = uniqueinds(new_data[i], new_data[i - 1])
            U, S, V = svd(core, right_inds...; cutoff=cutoff, maxdim=maxdim)

            new_data[i] = U
            new_data[i - 1] = S * V

            # Update bond
            bond = commonind(U, S)
            # Bond is between i and i-1.
            # i-1 is the left tensor.
            idx_left = i-1
            if isodd(idx_left)
                m = (idx_left+1)÷2
                new_bonds_copy[m] = bond
            else
                c = idx_left÷2
                if c <= length(new_bonds_main)
                    new_bonds_main[c] = bond
                end
            end
        end
        oc = 1
    else
        throw(
            ArgumentError(
                "zip_to_compress_mpo: Unknown direction '$direction'. Must be 'up' or 'down'.",
            ),
        )
    end

    return PairedSiteMPO(
        new_data, mpo.sites_main, mpo.sites_copy, new_bonds_main, new_bonds_copy
    ),
    oc
end

"""
    build_dt_mpo(n::Int, ωr::Real, sites_main, sites_copy; cutoff=1e-14, maxdim=1000)

Build the full Damping Transform MPO for n qubits with damping parameter ωr.

The DT MPO is composed of two parts:
1. Part 1: Blocks with control on main sites (k = 1 to n)
2. Part 2: Blocks with control on copy sites (k = n to 2)

Uses zip-to-combine and zip-to-compress algorithms following the paper.

Returns a PairedSiteMPO acting on 2n sites (n main + n copy).
"""
function build_dt_mpo(
    n::Int,
    ωr::Real,
    sites_main::Vector{I},
    sites_copy::Vector{I};
    cutoff=1e-14,
    maxdim=1000,
) where {I<:Index}
    n >= 1 || throw(ArgumentError("build_dt_mpo: n must be ≥ 1. Found n=$n"))
    length(sites_main) == n || throw(
        ArgumentError(
            "build_dt_mpo: sites_main must have $n elements. Got $(length(sites_main))"
        ),
    )
    length(sites_copy) == n || throw(
        ArgumentError(
            "build_dt_mpo: sites_copy must have $n elements. Got $(length(sites_copy))"
        ),
    )

    # Interleave sites: [main[1], copy[1], main[2], copy[2], ...]
    all_sites = Vector{I}(undef, 2n)
    for i in 1:n
        all_sites[2i - 1] = sites_main[i]
        all_sites[2i] = sites_copy[i]
    end

    if n == 1
        # For n=1, just Part 1 block for k=1 (damped Hadamard on main, identity on copy)
        return control_damping_mpo(1, 1, ωr, all_sites)
    end

    # =====================================================
    # Part 1: Build damping tensor train from l = 1 to k for k = 1 to n
    # =====================================================

    mpo_part1 = control_damping_mpo(n, 1, ωr, all_sites[1:2])
    oc = 0

    for k in 2:n
        # Extend mpo_part1 to include site k (initialized to Identity)
        # This is necessary because block_k acts on 1:k, while mpo_part1 currently acts on 1:k-1
        if length(mpo_part1.sites_main) < k
            # Get new sites
            s_main = all_sites[2 * k - 1]
            s_copy = all_sites[2 * k]

            # Create new bonds
            b_main = Index(1, "bond-main-$(k-1)")
            b_copy = Index(1, "bond-copy-$k")

            # Update last tensor of mpo_part1 (copy[k-1]) to connect to new main[k]
            mpo_part1.data[end] *= ITensor(1.0, b_main)
            push!(mpo_part1.bonds_main, b_main)

            # Create new tensors for site k
            # main[k]: connects to b_main (left) and b_copy (right)
            T_main = delta(s_main, s_main') * ITensor(1.0, b_main) * ITensor(1.0, b_copy)

            # copy[k]: connects to b_copy (left)
            T_copy = delta(s_copy, s_copy') * ITensor(1.0, b_copy)

            # Update mpo_part1
            push!(mpo_part1.data, T_main)
            push!(mpo_part1.data, T_copy)
            push!(mpo_part1.sites_main, s_main)
            push!(mpo_part1.sites_copy, s_copy)
            push!(mpo_part1.bonds_copy, b_copy)
        end

        # Block for control on main[k], acts on sites 1:2k
        block_k = control_damping_mpo(n, k, ωr, all_sites[1:2k])

        # Combine with current MPO (zip-down since they share first sites)
        mpo_part1, oc, _ = zip_to_combine_mpos(mpo_part1, block_k, oc)

        # Compress
        mpo_part1, oc = zip_to_compress_mpo(
            mpo_part1, oc, "down"; cutoff=cutoff, maxdim=maxdim
        )
    end

    # =====================================================
    # Part 2: Build copy tensor train from ℓ = k+1 to n for k = 1 to n-1
    # =====================================================
    # Zip-combine blocks from k = 1 to n-1
    for k in 1:(n - 1)
        # Block for control on copy[k], acts on sites 2k-1:2n
        L = n - k
        block_k = control_damping_copy_mpo(n, k, ωr, all_sites[(2k - 1):2n])

        # Combine with current MPO (zip-up since they share last sites)
        mpo_part1, oc, _ = zip_to_combine_mpos(mpo_part1, block_k, oc)

        # Compress
        mpo_part1, oc = zip_to_compress_mpo(
            mpo_part1, oc, "up"; cutoff=cutoff, maxdim=maxdim
        )
    end
    return mpo_part1
end

"""
    build_dt_mpo(ψ::zTMPS, ωr::Real; cutoff=1e-14, maxdim=1000)

Build the DT MPO for a zTMPS signal.
"""
function build_dt_mpo(ψ::zTMPS, ωr::Real; cutoff=1e-14, maxdim=1000)
    n = length(ψ.sites_main)
    return build_dt_mpo(n, ωr, ψ.sites_main, ψ.sites_copy; cutoff=cutoff, maxdim=maxdim)
end

end # module DTTransform
