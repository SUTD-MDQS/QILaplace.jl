#src/transforms/zt_transformer.jl
# This module builds the discrete Laplace transform aka z-transform that acts on a given MPS to apply the damping transform circuit. 
# It used the dt_gates from DTGates and zt_gates from ZTGates to build the full zT MPO as a PairedSiteMPO that acts on zTMPS representations of signals. 

# The zT transformer has three parts: 
# - Part 1: control_damping_mpo blocks with control on main site k, targets on site 1...k-1
# - Part 2: control_damping_copy_mpo blocks with control on main site k, targets on site k+1...N
# Part 3: control_Hphase_ztmps_mpo blocks with control on copy site k, targets on site 1...k-1

# The algorithm composes these blocks using zip-up and zip-down compression. 

module ZTTransformer

using ITensors, Printf
using ..Mpo: PairedSiteMPO
using ..Mps: zTMPS
using ..DTGates: control_damping_mpo, control_damping_copy_mpo
using ..ZTGates: control_Hphase_ztmps_mpo
using ..DTTransform: zip_to_combine_mpos, zip_to_compress_mpo

export build_zt_mpo

function build_zt_mpo(n::Int, ωr::Real, sites_main::Vector{I}, sites_copy::Vector{I};
                            cutoff=1e-14, maxdim=1000) where {I<:Index}
    n >= 1 || throw(ArgumentError("build_zt_mpo: n must be ≥ 1. Found n=$n"))
    length(sites_main) == n || throw(ArgumentError("build_zt_mpo: sites_main must have $n elements. Got $(length(sites_main))"))
    length(sites_copy) == n || throw(ArgumentError("build_zt_mpo: sites_copy must have $n elements. Got $(length(sites_copy))"))

    # Interleave sites: [main[1], copy[1], main[2], copy[2], ...]
    all_sites = Vector{I}(undef, 2n)
    for i in 1:n
        all_sites[2i - 1] = sites_main[i]
        all_sites[2i] = sites_copy[i]
    end

    if n == 1
        # For n=1, just part 1 block for k=1 (damped Hadamard on main, identity on copy) and part 3 block for k=1 (identity on main, Hadamard on copy)
        mpo_part1 = control_damping_mpo(1, 1, ωr, all_sites)
        mpo_part3 = control_Hphase_ztmps_mpo(1, all_sites)
        
        full_mpo, _, _ = zip_to_combine_mpos(mpo_part1, mpo_part3, 0)
        return full_mpo
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
            s_main = all_sites[2*k - 1]
            s_copy = all_sites[2*k]
            
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
        mpo_part1, oc = zip_to_compress_mpo(mpo_part1, oc, "down"; cutoff=cutoff, maxdim=maxdim)
    end

    # =====================================================
    # Part 2: Build copy tensor train from ℓ = k+1 to n for k = 1 down to n-1
    # =====================================================
    # Zip-combine blocks from k = 1 to n-1
    for k in 1:(n-1)
        # Block for control on main[k], acts on sites 2k+1:2n
        L = n - k
        block_k = control_damping_copy_mpo(n, k, ωr, all_sites[2k-1:2n])
        
        # Combine with current MPO (zip-up since they share last sites)
        mpo_part1, oc = zip_to_combine_mpos(mpo_part1, block_k, oc)
        
        # Compress
        mpo_part1, oc = zip_to_compress_mpo(mpo_part1, oc, "up"; cutoff=cutoff, maxdim=maxdim)
    end

    # =====================================================
    # Part 3: Build qft tensor train from k = 1 to n
    # =====================================================
    # Zip-combine blocks from k = 1 to n
    for k in 1:n
        # Block for control on copy[k], acts on sites 1:2k
        block_qft_k = control_Hphase_ztmps_mpo(k, all_sites[1:2k])

        # Combine with current MPO (zip-down since they share first sites)
        mpo_part1, oc = zip_to_combine_mpos(mpo_part1, block_qft_k, oc)

        # Compress
        mpo_part1, oc = zip_to_compress_mpo(mpo_part1, oc, "down"; cutoff=cutoff, maxdim=maxdim)
    end
    return mpo_part1
end

build_zt_mpo(ψ::zTMPS, ωr::Real; cutoff=1e-14, maxdim=1000) = build_zt_mpo(length(ψ.sites_main), ωr, ψ.sites_main, ψ.sites_copy; cutoff=cutoff, maxdim=maxdim)

end # module ZTTransformer
