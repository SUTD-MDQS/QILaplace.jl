using Test, QILaplace
using ITensors, LinearAlgebra, Printf


# ------------------------ helper functions ------------------------

################# SIGNAL CONVERTERS TEST HELPERS #################

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
    n_target >= n_window || throw(ArgumentError("embed_mpo: target length must be >= MPO length"))

    start = findfirst(==(W.sites[1]), target_sites)
    start === nothing && throw(ArgumentError("embed_mpo: window not found in target sites"))
    (start + n_window - 1) <= n_target || throw(ArgumentError("embed_mpo: window exceeds target range"))
    target_sites[start:start + n_window - 1] == W.sites || throw(ArgumentError("embed_mpo: target window must match MPO sites"))

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

################### QFT GATES TEST HELPERS ###################

_int_to_bit(b::Int, n::Int) = reverse(digits(b, base=2, pad=n))

