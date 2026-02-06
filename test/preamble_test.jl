using Test, QILaplace, ITensors
using Random, LinearAlgebra, Printf
using FFTW

# ------------------------ helper functions ------------------------

################# SIGNAL CONVERTERS TEST HELPERS #################

"""
Convert an integer to a bit vector of length n.
order=:msb -> [bit_n, ..., bit_1] (MSB first)
order=:lsb -> [bit_1, ..., bit_n] (LSB first)
"""
function int_to_bits(val::Int, n::Int; order::Symbol=:msb)
    if order == :msb
        return reverse(digits(val; base=2, pad=n))
    elseif order == :lsb
        return digits(val; base=2, pad=n)
    else
        throw(ArgumentError("order must be :msb or :lsb"))
    end
end

"""Convert a bit vector back to an integer."""
function bits_to_int(bits::AbstractVector{<:Integer}; order::Symbol=:msb)
    val = 0
    if order == :msb
        for b in bits
            val = (val << 1) | (b & 1)
        end
    elseif order == :lsb
        for (i, b) in enumerate(bits)
            val |= (b & 1) << (i - 1)
        end
    else
        throw(ArgumentError("order must be :msb or :lsb"))
    end
    return val
end

"""
Convert integer to interleaved main/copy bits for 2n-site representation.
Matches the interleaved ordering [main[1], copy[1], main[2], copy[2], ...].
"""
function int_to_paired_bits(val::Int, n::Int; order::Symbol=:msb)
    bits = int_to_bits(val, n; order=order)
    paired = Vector{Int}(undef, 2n)
    for i in 1:n
        paired[2i - 1] = bits[i]
        paired[2i] = bits[i]
    end
    return paired
end

function contract_chain(data::Vector{ITensor})
    T = ITensor(1)
    for A in data
        T *= A
    end
    return T
end

################### APPLY MPO TEST HELPERS ###################

function to_dense_mps(ψ::SignalMPS)
    n = length(ψ)
    result = ψ.data[1]
    for i in 2:n
        result = result * ψ.data[i]
    end
    return result
end

function to_dense_mpo(W::SingleSiteMPO)
    n = length(W)
    result = W.data[1]
    for i in 2:n
        result = result * W.data[i]
    end
    return result
end

function apply_dense(W::SingleSiteMPO, ψ::SignalMPS)
    W_dense = to_dense_mpo(W)
    ψ_dense = to_dense_mps(ψ)

    # Contract: MPO primed indices with MPS unprimed indices
    # First, replace primed MPO indices with temp indices
    result_mpo = copy(W_dense)
    result_mps = copy(ψ_dense)
    for site in ψ.sites
        site_tmp = sim(site)
        result_mpo = replaceind(result_mpo, site', site_tmp)
        result_mps = replaceind(result_mps, site, site_tmp)
    end

    # Contract
    result = result_mpo * result_mps

    return result
end

function dense_compose_mpos(W1::SingleSiteMPO, W2::SingleSiteMPO)
    dense1 = to_dense_mpo(W1)
    dense2 = to_dense_mpo(W2)

    # Step 1: lift W1 input (primed) indices to double-prime level
    for s in W1.sites
        dense1 = replaceind(dense1, prime(s), prime(s, 2))
    end

    # Step 2: lift W1 output indices to prime level so they match W2 inputs
    for s in W1.sites
        dense1 = replaceind(dense1, s, prime(s))
    end

    dense = dense1 * dense2

    # Step 3: restore double-prime indices back to single-prime to match MPO format
    for s in W1.sites
        dense = replaceind(dense, prime(s, 2), prime(s))
    end

    return dense
end

# Embed a SingleSiteMPO into a larger set of target sites
function embed_mpo(W::SingleSiteMPO, target_sites::Vector{<:Index})
    n_target = length(target_sites)
    n_window = length(W)
    n_target >= n_window ||
        throw(ArgumentError("embed_mpo: target length must be >= MPO length"))

    start = findfirst(==(W.sites[1]), target_sites)
    start === nothing && throw(ArgumentError("embed_mpo: window not found in target sites"))
    (start + n_window - 1) <= n_target ||
        throw(ArgumentError("embed_mpo: window exceeds target range"))
    target_sites[start:(start + n_window - 1)] == W.sites ||
        throw(ArgumentError("embed_mpo: target window must match MPO sites"))

    IndexType = eltype(target_sites)
    new_bonds = IndexType[]
    if n_target > 1
        resize!(new_bonds, n_target - 1)
        for i in 1:length(new_bonds)
            if start <= i && i < start + n_window - 1
                new_bonds[i] = W.bonds[i - start + 1]
            else
                new_bonds[i] = Index(1, "embed-bond-$i")
            end
        end
    end

    new_data = Vector{ITensor}(undef, n_target)
    for i in 1:n_target
        s = target_sites[i]
        sp = prime(s)
        if i < start || i > start + n_window - 1
            T = delta(sp, s)
            if i > 1
                T *= ITensor([1.0], new_bonds[i - 1])
            end
            if i <= length(new_bonds)
                T *= ITensor([1.0], new_bonds[i])
            end
            new_data[i] = T
        else
            w_idx = i - start + 1
            T = copy(W.data[w_idx])
            if w_idx == 1 && i > 1
                T *= ITensor([1.0], new_bonds[i - 1])
            end
            if w_idx == n_window && i <= length(new_bonds)
                T *= ITensor([1.0], new_bonds[i])
            end
            new_data[i] = T
        end
    end

    return SingleSiteMPO(new_data, target_sites, new_bonds)
end

############################### VECTOR EXTRACTION HELPERS ###############################

"""Create a computational-basis ket vector for the given bits."""
function basis_state_vector(bits::Vector{Int})
    n = length(bits)
    N = 2^n
    vec = zeros(N)
    index = sum(bits[i] * 2^(n - i) for i in 1:n) + 1
    vec[index] = 1.0
    return vec
end

"""Create a computational-basis ket vector for integer index i (0-indexed)."""
function basis_state_vector(i::Int, n::Int)
    N = 2^n
    @assert 0 <= i < N "Basis state index out of range"
    vec = zeros(Float64, N)
    vec[i + 1] = 1.0
    return vec
end

"""Extract state vector from SignalMPS (MSB-first bit ordering)."""
function mps_to_vector(ψ::SignalMPS)
    n = length(ψ.sites)
    N = 2^n

    T = ψ.data[1]
    for i in 2:n
        T *= ψ.data[i]
    end

    # Extract as array with proper ordering (MSB-first: reverse sites)
    arr = Array(T, reverse(ψ.sites)...)
    vec = reshape(arr, N)
    return vec
end

"""Extract state vector from zTMPS (flattens 2n-site representation)."""
function mps_to_vector(ψ::zTMPS)
    ψ2n = _as_signal_2n(ψ)
    T = prod(ψ2n.data)
    # Note: sites in ψ2n are interleaved [main[1], copy[1], main[2], copy[2], ...]
    # Extract as flat vector (LSB-first: site[1] is fastest-varying)
    A = Array(T, ψ2n.sites...)
    return ComplexF64.(vec(A))
end
