# src/Mpo.jl
# This module contains the structures and functions related to Matrix Product Operators (MPO) that will be used in the other modules involving constructing the transform MPO and applying them to MPS.
module Mpo

# Extend the generic interface functions for MPS/MPO index updates
import ..Mps: update_site!, update_bond!
using ITensors, Random, Printf

export SingleSiteMPO, PairedSiteMPO,
       check_singlesitempo, check_pairedsitempo

##################################### MPO TYPES #########################################

abstract type AbstractMPO end

# SingleSiteMPO can contain the Quantum Fourier Transform as an MPO in a n-qubit system
mutable struct SingleSiteMPO{I<:Index} <: AbstractMPO
    data::Vector{ITensor}
    sites::Vector{I}
    bonds::Vector{I}
    function SingleSiteMPO(data::Vector{ITensor}, sites::Vector{I}, bonds::Vector{I}) where I<:Index
        try 
            check_singlesitempo(data, sites, bonds)
        catch e
            rethrow(e)
        end
        new{I}(data, sites, bonds)
    end
end

# PairedSiteMPO can contain the z-transform and Exponential Decay Transform as an MPO in a 2n-qubit system
mutable struct PairedSiteMPO{I<:Index} <:AbstractMPO
    data::Vector{ITensor}
    sites_main::Vector{I}
    sites_copy::Vector{I}
    bonds_main::Vector{I}
    bonds_copy::Vector{I}
    function PairedSiteMPO(data::Vector{ITensor}, sites_main::Vector{I}, sites_copy::Vector{I}, bonds_main::Vector{I}, bonds_copy::Vector{I}) where I<:Index
        try
            check_pairedsitempo(data, sites_main, sites_copy, bonds_main, bonds_copy)
        catch e
            rethrow(e)
        end
        new{I}(data, sites_main, sites_copy, bonds_main, bonds_copy)
    end
end

##################################### MPO IDENTITY CONSTRUCTORS #########################################
function SingleSiteMPO(n::Int)
    sites = [Index(2, @sprintf("site-%d", i)) for i in 1:n]
    bonds = [Index(1, @sprintf("bond-%d", i)) for i in 1:(n - 1)]
    if n == 1
        return SingleSiteMPO([delta(sites[1], sites[1]')], sites, eltype(sites)[])
    end
    data = Vector{ITensor}(undef, n)
    for i in eachindex(data)
        if i == 1
            data[i] = delta(sites[i], sites[i]') * ITensor([1], bonds[i])
        elseif i == n
            data[i] = delta(sites[i], sites[i]') * ITensor([1], bonds[i - 1])
        else
            data[i] = delta(sites[i], sites[i]') * ITensor([1], bonds[i - 1]) * ITensor([1], bonds[i])
        end
    end
    return SingleSiteMPO(data, sites, bonds)
end

function PairedSiteMPO(n::Int)
    sites_main = [Index(2, @sprintf("site-main-%d", i)) for i in 1:n]
    sites_copy = [Index(2, @sprintf("site-main-%d", i)) for i in 1:n]
    bonds_main = [Index(1, @sprintf("bond-main-%d", i)) for i in 1:(n - 1)]
    bonds_copy = [Index(1, @sprintf("bond-copy-%d", i)) for i in 1:n]

    data = Vector{ITensor}(undef, 2n)

    for i in eachindex(data)
        if isodd(i) # odd are main sites
            m = div(i + 1, 2)
            if m == 1
                # main1: only intra bond copy[1]
                tb = ITensor([1], bonds_copy[m])
                data[i] = delta(sites_main[m], sites_main[m]') * tb
            elseif m == n
                # last main: left main[n-1], right copy[n]
                tbL = ITensor([1], bonds_main[m - 1])
                tbR = ITensor([1], bonds_copy[m])
                data[i] = delta(sites_main[m], sites_main[m]') * tbL * tbR
            else
                tbL = ITensor([1], bonds_main[m - 1])
                tbR = ITensor([1], bonds_copy[m])
                data[i] = delta(sites_main[m], sites_main[m]') * tbL * tbR
            end
        else # even are copy sites
            c = div(i, 2)
            if c == 1
                # first copy: left copy[1], right main[1] (if present)
                tbL = ITensor([1], bonds_copy[c])
                if n > 1
                    tbR = ITensor([1], bonds_main[c])
                    data[i] = delta(sites_copy[c], sites_copy[c]') * tbL * tbR
                else
                    data[i] = delta(sites_copy[c], sites_copy[c]') * tbL
                end
            elseif c == n
                # last copy: only intra bond copy[n]
                tb = ITensor([1], bonds_copy[c])
                data[i] = delta(sites_copy[c], sites_copy[c]') * tb
            else
                tbL = ITensor([1], bonds_main[c])
                tbR = ITensor([1], bonds_copy[c])
                data[i] = delta(sites_copy[c], sites_copy[c]') * tbL * tbR
            end
        end
    end
    return PairedSiteMPO(data, sites_main, sites_copy, bonds_main, bonds_copy)
end

####################################### MPO VALIDATION #########################################
function check_singlesitempo(data::Vector{ITensor}, sites::Vector{<:Index}, bonds::Vector{<:Index})
    n = length(sites)
    length(data) == n || throw(ArgumentError("SingleSiteMPO: Data length must equal number of sites. Got $n sites for $(length(data)) tensors."))
    length(bonds) == n - 1 || throw(ArgumentError("SingleSiteMPO: Number of bonds must be one less than number of sites. Got $(length(bonds)) bonds for $n sites."))
    
    if n == 1
        length(inds(data[1])) == 2 || throw(ArgumentError("SingleSiteMPO: Data tensor at site 1 must have exactly 2 indices. Found $(length(inds(data[1])))"))
        (sites[1] in inds(data[1])) || throw(ArgumentError("SingleSiteMPO: Site index missing in single-site tensor"))
        return nothing
    end

    # Edge ranks = 3, bulk ranks = 4
    (length(inds(data[1])) == 3) || throw(ArgumentError("SingleSiteMPO: Edge tesors must have exactly 3 indices. Found $(length(inds(data[1]))) at site 1"))
    (length(inds(data[n])) == 3) || throw(ArgumentError("SingleSiteMPO: Edge tesors must have exactly 3 indices. Found $(length(inds(data[n]))) at site $n"))
    for i in 2:(n - 1)
        (length(inds(data[i])) == 4) || throw(ArgumentError("SingleSiteMPO: Bulk tesors must have exactly 4 indices. Found $(length(inds(data[i]))) at site $i"))
    end

    # Site and its prime presence and uniqueness
    for i in 1:n
        (sites[i] in inds(data[i]) && sites[i]' in inds(data[i])) || throw(ArgumentError("SingleSiteMPO: Site index or its prime missing in tensor at site $i. Got indices: $(inds(data[i]))"))
    end
    length(unique(sites)) == n || throw(ArgumentError("SingleSiteMPO: Site indices must be unique"))

    # Bond presence and uniqueness
    for i in 1:(n - 1)
        (bonds[i] in inds(data[i])) || throw(ArgumentError("SingleSiteMPO: Bond index missing in tensor at site $i"))
        (bonds[i] === commonind(data[i], data[i + 1])) || throw(ArgumentError("SingleSiteMPO: Bond mismatch between core at $i and $(i + 1)"))
    end
    length(unique(bonds)) == n - 1 || throw(ArgumentError("SingleSiteMPO: Bond indices must be unique"))
    return nothing
end
# Convenience checker 
check_singlesitempo(W::SingleSiteMPO) = check_singlesitempo(W.data, W.sites, W.bonds)

# Length and iteration protocol for SingleSiteMPO
Base.length(W::SingleSiteMPO) = length(W.data)
Base.iterate(W::SingleSiteMPO) = iterate(W.data)
Base.iterate(W::SingleSiteMPO, s) = iterate(W.data, s)

function check_pairedsitempo(data::Vector{ITensor}, sites_main::Vector{<:Index}, sites_copy::Vector{<:Index}, bonds_main::Vector{<:Index}, bonds_copy::Vector{<:Index})
    n = length(sites_main)
    length(data) == 2n || throw(ArgumentError("PairedSiteMPO: Data length must equal twice the number of sites. Got $n sites for $(length(data)) tensors."))
    length(bonds_main) == n - 1 || throw(ArgumentError("PairedSiteMPO: Number of main bonds must be n-1. Got $(length(bonds_main)) bonds for $n sites."))
    length(bonds_copy) == n || throw(ArgumentError("PairedSiteMPO: Number of copy bonds must be n. Got $(length(bonds_copy)) bonds for $n sites."))

    # Check if the main and copy sites/bonds are unique and dont have overlap between them
    commoninds = intersect(sites_main, sites_copy)
    length(commoninds) == 0 || throw(ArgumentError("PairedSiteMPO: sites_main and sites_copy must be disjoint sets. Found common indices: $commoninds"))
    commoninds = intersect(bonds_main, bonds_copy)
    length(commoninds) == 0 || throw(ArgumentError("PairedSiteMPO: bonds_main and bonds_copy must be disjoint sets. Found common indices: $commoninds"))
    
    if n == 1
        length(inds(data[1])) == 3 || throw(ArgumentError("PairedSiteMPO: Main data tensor at site 1 must have exactly 3 indices. Found $(length(inds(data[1])))"))
        (sites_main[1] in inds(data[1]) && sites_main[1]' in inds(data[1])) || throw(ArgumentError("PairedSiteMPO: Main site index or its prime missing in single-site tensor"))
        length(inds(data[2])) == 3 || throw(ArgumentError("PairedSiteMPO: Copy data tensor at site 1 must have exactly 3 indices. Found $(length(inds(data[2])))"))
        (sites_copy[1] in inds(data[2]) && sites_copy[1]' in inds(data[2])) || throw(ArgumentError("PairedSiteMPO: Copy site index or its prime missing in single-site tensor"))
        commonind(data[1], data[2]) === bonds_copy[1] || throw(ArgumentError("PairedSiteMPO: Copy bond mismatch between main and copy tensor at site 1. Got $(commonind(data[1], data[2])) but expected $(bonds_copy[1])"))
        return nothing
    end

    # Edge ranks = 3, bulk ranks = 4
    for i in 1:2n
        if i == 1 || i == 2n
            (length(inds(data[i])) == 3) || throw(ArgumentError("PairedSiteMPO: Edge tesors must have exactly 3 indices. Found $(length(inds(data[i]))) at site $i"))
        else
            (length(inds(data[i])) == 4) || throw(ArgumentError("PairedSiteMPO: Bulk tesors must have exactly 4 indices. Found $(length(inds(data[i]))) at site $i"))
        end
    end

    # Site and its prime presence and uniqueness
    for i in 1:2n
        if isodd(i) # odd are main sites
            site_idx = div(i + 1, 2)
            (sites_main[site_idx] in inds(data[i]) && sites_main[site_idx]' in inds(data[i])) || throw(ArgumentError("PairedSiteMPO: Site index or its prime missing in tensor at site $i. Got indices: $(inds(data[i]))"))
        else # even are copy sites
            site_idx = div(i, 2)
            (sites_copy[site_idx] in inds(data[i]) && sites_copy[site_idx]' in inds(data[i])) || throw(ArgumentError("PairedSiteMPO: Site index or its prime missing in tensor at site $i. Got indices: $(inds(data[i]))"))
        end
    end
    length(unique(sites_main)) == n || throw(ArgumentError("PairedSiteMPO: Site indices must be unique in sites_main"))
    length(unique(sites_copy)) == n || throw(ArgumentError("PairedSiteMPO: Site indices must be unique in sites_copy"))


    # Bond presence and uniqueness
    for i in 1:(2n - 1)
        if isodd(i) # odd bonds are copy bonds between main(i) and copy(i)
            bond_idx = div(i + 1, 2)
            (bonds_copy[bond_idx] in inds(data[i])) || throw(ArgumentError("PairedSiteMPO: Copy bond index missing in tensor at site $i"))
            (bonds_copy[bond_idx] === commonind(data[i], data[i + 1])) || throw(ArgumentError("PairedSiteMPO: Copy bond mismatch between core at $i and $(i + 1)"))
        else # even bonds are main bonds between copy(i) and main(i+1)
            bond_idx = div(i, 2)
            (bonds_main[bond_idx] in inds(data[i])) || throw(ArgumentError("PairedSiteMPO: Main bond index missing in tensor at site $i"))
            (bonds_main[bond_idx] === commonind(data[i], data[i + 1])) || throw(ArgumentError("PairedSiteMPO: Main bond mismatch between core at $i and $(i + 1)"))
        end
    end
    length(unique(bonds_main)) == n - 1 || throw(ArgumentError("PairedSiteMPO: Bond indices must be unique in bonds_main"))
    length(unique(bonds_copy)) == n || throw(ArgumentError("PairedSiteMPO: Bond indices must be unique in bonds_copy"))
    return nothing
end
# Convenience checker
check_pairedsitempo(W::PairedSiteMPO) = check_pairedsitempo(W.data, W.sites_main, W.sites_copy, W.bonds_main, W.bonds_copy)


###################################### DISPLAY UTILS #########################################
function Base.show(io::IO, W::SingleSiteMPO)
    n = length(W.sites)
    println(io, "SingleSiteMPO with $n sites:")
    for i in 1:n
        t = W.data[i]
        idxs = collect(inds(t))
        parts = String[]
        for I in idxs
            d = dim(I)
            tg = try s = string(tags(I)); isempty(s) ? nothing : s catch; nothing end
            push!(parts, tg === nothing ? "dim=$d" : "dim=$d, tags=$tg")
        end
        println(io, "  Site $i: ", join(parts, " | "))
    end
end

function Base.show(io::IO, W::PairedSiteMPO)
    n = length(W.sites_main)
    println(io, "PairedSiteMPO with $n sites:")
    for i in 1:n
        # Main site tensor (odd index in data)
        t_main = W.data[2i - 1]
        idxs_main = collect(inds(t_main))
        parts_main = String[]
        for I in idxs_main
            d = dim(I)
            tg = try s = string(tags(I)); isempty(s) ? nothing : s catch; nothing end
            push!(parts_main, tg === nothing ? "dim=$d" : "dim=$d, tags=$tg")
        end
        println(io, "  Site $i (main): ", join(parts_main, " | "))
        
        # Copy site tensor (even index in data)
        t_copy = W.data[2i]
        idxs_copy = collect(inds(t_copy))
        parts_copy = String[]
        for I in idxs_copy
            d = dim(I)
            tg = try s = string(tags(I)); isempty(s) ? nothing : s catch; nothing end
            push!(parts_copy, tg === nothing ? "dim=$d" : "dim=$d, tags=$tg")
        end
        println(io, "  Site $i (copy): ", join(parts_copy, " | "))
    end
end

###################################### MPO UPDATE HELPERS #########################################

# Update a site by specifying old and new Index values for SingleSiteMPO
function update_site!(W::SingleSiteMPO, old_site_index::Index, new_site_index::Index)
    # find site position
    i = findfirst(x -> x == old_site_index, W.sites)
    if i === nothing
        # fallback: search data tensors for either old or its prime
        i = findfirst(k -> (old_site_index in inds(W.data[k]) || prime(old_site_index) in inds(W.data[k])), 1:length(W.data))
        i === nothing && throw(ArgumentError("update_site!: old site index not found in SingleSiteMPO"))
    end

    E = typeof(W.sites[i])
    (new_site_index isa E) || throw(ArgumentError("update_site!: new site index must be of type $E"))
    dim(old_site_index) == dim(new_site_index) || throw(ArgumentError("update_site!: Site index dimension mismatch at site $i"))

    # replace unprimed and primed occurrences
    replaceinds!(W.data[i], old_site_index => new_site_index)
    replaceinds!(W.data[i], prime(old_site_index) => prime(new_site_index))
    W.sites[i] = new_site_index
    return W
end

# Update a bond by specifying old and new Index values for SingleSiteMPO
function update_bond!(W::SingleSiteMPO, old_bond_index::Index, new_bond_index::Index)
    j = findfirst(x -> x == old_bond_index, W.bonds)
    j === nothing && throw(ArgumentError("update_bond!: old bond index not found in SingleSiteMPO"))
    B = typeof(W.bonds[j])
    (new_bond_index isa B) || throw(ArgumentError("update_bond!: new bond index must be of type $B"))
    dim(old_bond_index) == dim(new_bond_index) || throw(ArgumentError("update_bond!: Bond index dimension mismatch at bond $j"))

    replaceinds!(W.data[j], old_bond_index => new_bond_index)
    replaceinds!(W.data[j + 1], old_bond_index => new_bond_index)
    W.bonds[j] = new_bond_index
    return W
end

# Update a site by specifying old and new Index values for PairedSiteMPO
function update_site!(W::PairedSiteMPO, old_site_index::Index, new_site_index::Index)
    m = findfirst(x -> x == old_site_index, W.sites_main)
    if m !== nothing
        E = typeof(W.sites_main[m])
        (new_site_index isa E) || throw(ArgumentError("update_site!: new main site index must be of type $E"))
        dim(old_site_index) == dim(new_site_index) || throw(ArgumentError("update_site!: Main site index dimension mismatch at site $m"))
        replaceinds!(W.data[2m - 1], old_site_index => new_site_index)
        replaceinds!(W.data[2m - 1], prime(old_site_index) => prime(new_site_index))
        W.sites_main[m] = new_site_index
        return W
    end

    c = findfirst(x -> x == old_site_index, W.sites_copy)
    if c !== nothing
        E = typeof(W.sites_copy[c])
        (new_site_index isa E) || throw(ArgumentError("update_site!: new copy site index must be of type $E"))
        dim(old_site_index) == dim(new_site_index) || throw(ArgumentError("update_site!: Copy site index dimension mismatch at site $c"))
        replaceinds!(W.data[2c], old_site_index => new_site_index)
        replaceinds!(W.data[2c], prime(old_site_index) => prime(new_site_index))
        W.sites_copy[c] = new_site_index
        return W
    end

    throw(ArgumentError("update_site!: old site index not found in PairedSiteMPO"))
end

# Update a bond by specifying old and new Index values for PairedSiteMPO
function update_bond!(W::PairedSiteMPO, old_bond_index::Index, new_bond_index::Index)
    # main bonds: connect copy(k) (data[2k]) and main(k+1) (data[2k+1])
    k = findfirst(x -> x == old_bond_index, W.bonds_main)
    if k !== nothing
        B = typeof(W.bonds_main[k])
        (new_bond_index isa B) || throw(ArgumentError("update_bond!: new main bond index must be of type $B"))
        dim(old_bond_index) == dim(new_bond_index) || throw(ArgumentError("update_bond!: Main bond dimension mismatch at bond $k"))
        replaceinds!(W.data[2k], old_bond_index => new_bond_index)
        replaceinds!(W.data[2k + 1], old_bond_index => new_bond_index)
        W.bonds_main[k] = new_bond_index
        return W
    end

    # copy bonds: connect main(j) (data[2j-1]) and copy(j) (data[2j])
    j = findfirst(x -> x == old_bond_index, W.bonds_copy)
    if j !== nothing
        B = typeof(W.bonds_copy[j])
        (new_bond_index isa B) || throw(ArgumentError("update_bond!: new copy bond index must be of type $B"))
        dim(old_bond_index) == dim(new_bond_index) || throw(ArgumentError("update_bond!: Copy bond dimension mismatch at bond $j"))
        replaceinds!(W.data[2j - 1], old_bond_index => new_bond_index)
        replaceinds!(W.data[2j], old_bond_index => new_bond_index)
        W.bonds_copy[j] = new_bond_index
        return W
    end

    throw(ArgumentError("update_bond!: old bond index not found in PairedSiteMPO"))
end

end # module Mpo
