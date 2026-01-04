# src/signals/SignalConverters.jl
# This module provides functions to convert signals (1D arrays) into MPS representations using SVD and RSVD methods.

module SignalConverters

using ITensors, Random, Printf
using LinearAlgebra: norm
using ..RSVD: rsvd

using ..Mps: AbstractMPS, SignalMPS, zTMPS,
             PairCore,
             nsite, siteindices, bondindices

export signal_mps, signal_ztmps

################################## SIGNAL VEC TO MPS #########################################
# Convert an input 1D array signal into ITensor MPS with the binary encoding of the indices
function _array_to_tensor(x::AbstractVector; sites = undef)
    N = length(x)
    n = round(Int, log2(N)) # no. of qubits that can encodes signal

    ITensors.disable_warn_order()
    try
        # If signal is not a power of 2, fill with 0s upto length 2^n with a warning
        if N < 2^n
            @warn "_array_to_tensor: Input signal length $N is not a power of 2. Filling with zeros upto length $(2^n). We recommend providing signals of length power-of-2 for best performance."
            x_filled = zeros(2^n)
            x_filled[1:N] .= x
            x = x_filled
        end
        @assert length(x) == 2^n || "_array_to_tensor: Length of signal vector must be a power of 2"
        if sites === undef
            sites = [Index(2; tags=@sprintf("site-%d", i)) for i in 1:n]
        end
        @assert length(sites) == n || "_array_to_tensor: Number of provided sites must match log2(length(x))"
        normalisation_const = norm(x)
        x /= normalisation_const

        x_tensor = reshape(x, ntuple(_ -> 2, n)...)
        ψ = ITensor(x_tensor, reverse(sites)...)
        ψ = permute(ψ, sites...) # first index in inds(ψ) corresponds to site-1 and so on...
        return ψ, normalisation_const
    finally
        ITensors.reset_warn_order()
    end
end

# Convert an ITensor to MPS with given tol. TODO: test this with a vector of [1 2 3 ... 64] and see if it transforms as expected
function _tensor_to_mps_svd(ψ::ITensor;
                      cutoff::Real = 1e-15, maxdim::Int = typemax(Int))
    ITensors.disable_warn_order()
    try
        N = length(inds(ψ))
        I = eltype(inds(ψ))
        
        N ≥ 1 || throw(ArgumentError("_tensor_to_mps_svd: Need at least one site in the tensor to convert to MPS. Found $N sites."))

        all(dims(inds(ψ)) .== 2) || throw(ArgumentError("_tensor_to_mps_svd: All indices of the input tensor must be of dimension 2 (qubit indices). Found dims: $(dims(inds(ψ)))"))

        # Trivial case: single site MPS
        if N == 1
            return SignalMPS([ψ], collect(inds(ψ)), Vector{I}(undef, 0))
        end
        
        data  = Vector{ITensor}(undef, N)
        sites = collect(inds(ψ))
        bonds = Vector{I}(undef, N - 1)
        current_ψ = ψ

        for i in 1:(N - 1)
            if i == 1
                left_inds = sites[i]
            else
                left_inds = (sites[i], bonds[i - 1])
            end

            U, S, V = svd(current_ψ, left_inds; cutoff = cutoff, maxdim = maxdim)
            old_bond = commonind(U, S)
            @assert !isnothing(old_bond) || error("No shared bond index between U and S at site $i")
            new_bond = Index(dim(old_bond); tags=@sprintf("bond-%d", i))

            replaceind!(U, old_bond, new_bond)
            replaceind!(S, old_bond, new_bond)

            data[i] = U
            bonds[i] = new_bond
            current_ψ = S * V
            
            @assert commonind(data[i], current_ψ) == bonds[i] "Bond index mismatch at site $i during MPS conversion."
        end
        data[N] = current_ψ
        return SignalMPS(data, sites, bonds)
    finally
        ITensors.reset_warn_order()
    end
end

"""
    _tensor_to_mps_rsvd(ψ::ITensor; cutoff=1e-15, maxdim=typemax(Int), kwargs...)

Convert a tensor `ψ` to an MPS using a divide-and-conquer Randomized SVD (RSVD) approach.
Recursively splits the tensor into left and right halves, applying RSVD at the cut to reduce bond dimensions.
"""
function _tensor_to_mps_rsvd(ψ::ITensor; cutoff::Real = 1e-15, maxdim::Int = typemax(Int), kwargs...)
    ITensors.disable_warn_order()
    try
        N = length(inds(ψ))
        N ≥ 1 || throw(ArgumentError("_tensor_to_mps_rsvd: Need at least one site in the tensor to convert to MPS. Found $N sites."))
        all(dims(inds(ψ)) .== 2) || throw(ArgumentError("_tensor_to_mps_rsvd: All indices of the input tensor must be of dimension 2 (qubit indices). Found dims: $(dims(inds(ψ)))"))
        # Trivial case: single site MPS
        if N == 1
            return SignalMPS([ψ], collect(inds(ψ)), Vector{Index}(undef, 0))
        end

        sites = collect(inds(ψ))
        I = eltype(sites)

        # Forward cutoff/maxdim via kwargs so rsvd gets them; they can be overridden by explicit kwargs
        kwargs = merge(values(kwargs), (; cutoff=cutoff, maxdim=maxdim))

        # Prepare containers
        data = Vector{ITensor}(undef, N)
        bonds = Vector{I}(undef, N - 1)

        # Attach trivial boundary links so we can consistently split and permute tensors
        left_bond_0 = Index(1, "dummy-bond-0")
        right_bond_N = Index(1, "dummy-bond-$N")
        T0 = ψ * ITensor(1.0, left_bond_0) * ITensor(1.0, right_bond_N)

        # Recursive compress function (inlined divide-and-conquer)
        function compress_tt!(T_chunk::ITensor, first_site_idx::Int, last_site_idx::Int, left_bond::Index, right_bond::Index)
            if first_site_idx == last_site_idx
                core = T_chunk
                if inds(core) != IndexSet(left_bond, sites[first_site_idx], right_bond)
                    core = permute(core, left_bond, sites[first_site_idx], right_bond)
                end
                data[first_site_idx] = core
                return
            end

            mid = (first_site_idx + last_site_idx - 1) ÷ 2
            left_inds = (left_bond, sites[first_site_idx:mid]...)

            # RSVD: Tseg ≈ U * S * V, where U has (lb, left sites..., uind)
            # and S*V has (uind, right sites..., rb).
            U, S, V = rsvd(T_chunk, left_inds...; kwargs...)
            T_left = U
            T_right = S * V

            # Record the link dimension between mid and mid+1.
            old_bond = commonind(T_left, T_right)
            @assert !isnothing(old_bond) || error("_tensor_to_mps_rsvd: missing internal bond between sites $mid and $(mid + 1)")
            new_bond = Index(dim(old_bond); tags=@sprintf("bond-%d", mid))
            bonds[mid] = new_bond

            replaceind!(T_left, old_bond, new_bond)
            replaceind!(T_right, old_bond, new_bond)

            compress_tt!(T_left, first_site_idx, mid, left_bond, new_bond)
            compress_tt!(T_right, mid + 1, last_site_idx, new_bond, right_bond)
            return
        end

        compress_tt!(T0, 1, N, left_bond_0, right_bond_N)

        # Strip boundary links so edge cores match the usual convention.
        data[1] = (data[1] * ITensor(1.0, left_bond_0))
        data[N] = (data[N] * ITensor(1.0, right_bond_N))

        return SignalMPS(data, sites, bonds)
    finally
        ITensors.reset_warn_order()
    end
end

# Simple dispatching wrapper that returns the same (data, bonds) pair expected by higher-level helpers
function _tensor_to_mps(ψ::ITensor; method::Symbol = :svd, kwargs...)
    method ∈ (:svd, :rsvd) || throw(ArgumentError("tensor_to_mps: unknown method $method. Use :svd or :rsvd."))
    if method == :svd
        return _tensor_to_mps_svd(ψ; kwargs...)
    else
        return _tensor_to_mps_rsvd(ψ; kwargs...)
    end
end

function signal_mps(x::AbstractVector{<:Number}; method::Symbol = :svd,
                    kwargs...)
    ψ, normalisation_const = _array_to_tensor(x)
    return _tensor_to_mps(ψ; method=method, kwargs...), normalisation_const
end

function signal_ztmps(x::AbstractVector{<:Number};
                         cutoff::Real = 1e-10, maxdim::Int = typemax(Int))
    # The SignalMPS to be copied
    ψ_signal, normalisation_const = signal_mps(x; cutoff=cutoff, maxdim=maxdim)
    n = nsite(ψ_signal)

    sites_main = sim.(ψ_signal.sites)
    sites_copy = [Index(dim(sites_main[i]); tags=@sprintf("site-copy-%d", i)) for i in 1:n]
    bonds_main = ψ_signal.bonds
    I = eltype(bonds_main)
    bonds_copy = Vector{I}(undef, n)

    paircores = Vector{PairCore}(undef, n)
    for i in 1:n
        # Build fused site tensors (project |σ⟩ -> |σ_main⟩ ⊗ |σ_copy⟩)
        T_core = ψ_signal.data[i] * delta(ψ_signal.sites[i], sites_main[i], sites_copy[i])

        left_inds = i == 1 ? (sites_main[i],) : (bonds_main[i - 1], sites_main[i])
        U, S, V = svd(T_core, left_inds...; cutoff=cutoff, maxdim=maxdim)
        core_main, core_copy = U, S * V

        old_bond = commonind(core_main, core_copy)
        @assert !isnothing(old_bond) || error("No shared bond index between main and copy cores at site $i")
        new_bond = Index(dim(old_bond); tags=@sprintf("bond-copy-%d", i))
        replaceind!(core_main, old_bond, new_bond)
        replaceind!(core_copy, old_bond, new_bond)

        paircores[i] = PairCore(core_main, core_copy, new_bond)
        bonds_copy[i] = new_bond
    end

    return zTMPS(paircores, bonds_main, bonds_copy, sites_main, sites_copy), normalisation_const
end

end # module SignalConverters
