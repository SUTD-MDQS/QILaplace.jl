# src/mps.jl
# This module contains the structures and functions related to Matrix Product States (MPS) that will be used in the main simulation code.
module MPS

using ITensors, Printf

export AbstractMPS, SignalMPS, QFTMPS, EDTMPS, zTMPS,
       PairCore,
       nsite, siteindices, bondindices,
       main_sites, copy_sites, main_bonds, copy_bonds,
       as_signal_2n, writeback_signal_2n!

##################################### MPS TYPES #########################################

abstract type AbstractMPS end

# Core structure for storing the two-leg zTMPS tensors
mutable struct PairCore
    Amain::ITensor
    Acopy::ITensor
    c::Index
    function PairCore(Amain::ITensor, Acopy::ITensor, c::Index)
        # sanity: Amain & Acopy must share `c`
        shared = commonind(inds(Amain), inds(Acopy))
        (shared === c) || throw(ArgumentError("Amain/Acopy must share bond c; got $(shared), expected $(c)"))
        new(Amain, Acopy, c)
    end
end

mutable struct SignalMPS{I<:Index} <: AbstractMPS
    data::Vector{ITensor{I}}
    sites::Vector{Index{I}}
    bonds::Vector{Index{I}}
    function SignalMPS(data::Vector{ITensor{I}}, sites::Vector{Index{I}}, bonds::Vector{Index{I}}) where I<:Index
        try 
            check_mps(data, sites, bonds)
        catch e
            throw(ArgumentError("Invalid SignalMPS construction: $(e.msg)"))
        end
        new{I}(data, sites, bonds)
    end
end

mutable struct QFTMPS{I<:Index} <: AbstractMPS
    ψ::SignalMPS{I}
end
mutable struct EDTMPS{I<:Index} <: AbstractMPS
    ψ::SignalMPS{I}
end

mutable struct zTMPS{I<:Index} <: AbstractMPS
    data::Vector{PairCore}
    bonds_main::Vector{Index{I}}
    bonds_copy::Vector{Index{I}}
    sites_main::Vector{Index{I}}
    sites_copy::Vector{Index{I}}

    function zTMPS(cores::Vector{PairCore},
                   bonds_main::Vector{Index{I}},
                   bonds_copy::Vector{Index{I}},
                   sites_main::Vector{Index{I}},
                   sites_copy::Vector{Index{I}}) where I<:Index
        try 
            check_ztmps(cores, bonds_main, bonds_copy, sites_main, sites_copy)
        catch e
            throw(ArgumentError("Invalid zTMPS construction: $(e.msg)"))
        end

        new(cores, bonds_main, bonds_copy, sites_main, sites_copy)
    end
end

##################################### MPS CONSTRUCTORS #########################################

function SignalMPS(sites::Vector{Index{I}}, bonds::Vector{Index{I}}) where I<:Index
    N = length(sites)
    N ≥ 1 || throw(ArgumentError("Need at least 1 site"))
    length(bonds) == N - 1 || throw(ArgumentError("Need N-1 bonds"))
    data = Vector{ITensor{I}}(undef, N)
    if N == 1
        data[1] = ITensor(sites[1])
    else
        data[1] = ITensor(sites[1], bonds[1])
        for i in 2:N-1
            data[i] = ITensor(bonds[i-1], sites[i], bonds[i])
        end
        data[N] = ITensor(bonds[N-1], sites[N])
    end
    SignalMPS(data, sites, bonds)
end

SignalMPS(N::Int, site_dim::Int) = begin
    sites = [Index(site_dim, @sprintf("site-%d", i)) for i in 1:N]
    bonds = [Index(1, @sprintf("bond-%d", i)) for i in 1:N-1]
    SignalMPS(sites, bonds)
end

function zTMPS(sites_main::Vector{I}, sites_copy::Vector{I}) where {I<:Index}
    N = length(sites_main)
    N == length(sites_copy) ||
        throw(ArgumentError("sites_main and sites_copy must have same length"))

    bonds_main = Vector{I}(undef, max(N - 1, 0))
    for i in eachindex(bonds_main)
        bonds_main[i] = Index(1, @sprintf("bond-main-%d", i))::I
    end

    bonds_copy = Vector{I}(undef, N)
    for i in 1:N
        bonds_copy[i] = Index(1, @sprintf("bond-copy-%d", i))::I
    end

    data = Vector{PairCore}(undef, N)

    for i in 1:N
        c = bonds_copy[i] # intra-site bond main(i) <-> copy(i)

        bL = i == 1 ? nothing : bonds_main[i - 1]
        bR = i == N ? nothing : bonds_main[i]

        Amain = bL === nothing ?
            ITensor(sites_main[i], c) :
            ITensor(bL, sites_main[i], c)
        Acopy = bR === nothing ?
            ITensor(c, sites_copy[i]) :
            ITensor(c, sites_copy[i], bR)

        data[i] = PairCore(Amain, Acopy, c)
    end

    return zTMPS(data, bonds_main, bonds_copy, sites_main, sites_copy)
end

zTMPS(N::Int, site_dim::Int) = begin
    sites_main = [Index(site_dim, @sprintf("site-main-%d", i)) for i in 1:N]
    sites_copy = [Index(site_dim, @sprintf("site-copy-%d", i)) for i in 1:N]
    zTMPS(sites_main, sites_copy)
end
##################################### MPS VALIDATION #########################################

function check_mps(data::Vector{ITensor}, sites::Vector{Index}, bonds::Vector{Index})
    N = length(data)
    length(sites) == N || throw(ArgumentError("Need $N site indices, got $(length(sites))"))
    length(bonds) == max(N - 1, 0) || throw(ArgumentError("Need $(N-1) bonds, got $(length(bonds))"))

    if N == 1
        @assert length(inds(data[1])) == 1 "Single-site MPS must have rank 1"
        (sites[1] in inds(data[1])) || throw(ArgumentError("Site index missing in single-site tensor"))
        return nothing
    end

    # Edge ranks = 2, bulk ranks = 3
    (length(inds(data[1])) == 2 && length(inds(data[N])) == 2) ||
        throw(ArgumentError("Edge tensors must have rank 2"))
    for i in 2:N-1
        length(inds(data[i])) == 3 || throw(ArgumentError("Bulk tensor at $i must have rank 3"))
    end

    # Site presence and uniqueness
    for i in 1:N
        (sites[i] in inds(data[i])) || throw(ArgumentError("Site index missing at site $i"))
    end
    length(unique(sites)) == N || throw(ArgumentError("Site indices must be unique"))

    # Bond wiring (open boundary)
    for i in 1:N-1
        shared = commonind(inds(data[i]), inds(data[i+1]))
        (shared === bonds[i]) || throw(ArgumentError("Bond mismatch between $i and $(i+1)"))
    end
    length(unique(bonds)) == length(bonds) || throw(ArgumentError("Bond indices must be unique"))
    return nothing
end

check_mps(ψ::SignalMPS) = check_mps(ψ.data, ψ.sites, ψ.bonds)

function check_ztmps(data::Vector{PairCore},
                     bonds_main::Vector{I},
                     bonds_copy::Vector{I},
                     sites_main::Vector{I},
                     sites_copy::Vector{I}) where {I<:Index}
    N = length(data)

    length(sites_main) == N ||
        throw(ArgumentError("zTMPS: sites_main must have length N"))
    length(sites_copy) == N ||
        throw(ArgumentError("zTMPS: sites_copy must have length N"))
    length(bonds_copy) == N ||
        throw(ArgumentError("zTMPS: bonds_copy (intra bonds) must have length N"))
    length(bonds_main) == max(N - 1, 0) ||
        throw(ArgumentError("zTMPS: bonds_main (inter bonds) must have length N-1"))

    # Optional uniqueness checks (can relax if needed)
    length(unique(sites_main)) == N ||
        throw(ArgumentError("zTMPS: sites_main must be unique"))
    length(unique(sites_copy)) == N ||
        throw(ArgumentError("zTMPS: sites_copy must be unique"))
    length(unique(bonds_copy)) == N ||
        throw(ArgumentError("zTMPS: bonds_copy (intra) must be unique"))
    if N > 1
        length(unique(bonds_main)) == length(bonds_main) ||
            throw(ArgumentError("zTMPS: bonds_main (inter) must be unique"))
    end

    for i in 1:N
        core = data[i]
        isA = inds(core.Amain)
        isB = inds(core.Acopy)

        # site indices present
        sites_main[i] in isA ||
            throw(ArgumentError("zTMPS: s_main[$i] missing in Amain[$i]"))
        sites_copy[i] in isB ||
            throw(ArgumentError("zTMPS: s_copy[$i] missing in Acopy[$i]"))

        # intra bond consistency
        bonds_copy[i] == core.c ||
            throw(ArgumentError("zTMPS: bonds_copy[$i] ≠ core.c at site $i"))
        core.c in isA ||
            throw(ArgumentError("zTMPS: intra bond c[$i] missing in Amain[$i]"))
        core.c in isB ||
            throw(ArgumentError("zTMPS: intra bond c[$i] missing in Acopy[$i]"))

        # left inter-site bond: copy(i-1) ↔ main(i)
        if i > 1
            bL = bonds_main[i - 1]
            bL in isA ||
                throw(ArgumentError("zTMPS: left inter bond b_main[$(i-1)] missing in Amain[$i]"))
            bL in inds(data[i - 1].Acopy) ||
                throw(ArgumentError("zTMPS: left inter bond b_main[$(i-1)] not in Acopy[$(i-1)]"))
        end

        # right inter-site bond: copy(i) ↔ main(i+1)
        if i < N
            bR = bonds_main[i]
            bR in isB ||
                throw(ArgumentError("zTMPS: right inter bond b_main[$i] missing in Acopy[$i]"))
        end
    end

    return nothing
end

check_ztmps(ψ::zTMPS) = check_ztmps(ψ.data, ψ.bonds_main, ψ.bonds_copy, ψ.sites_main, ψ.sites_copy)

##################################### LIGHTWEIGHT API HELPERS #########################################

nsite(ψ::SignalMPS) = length(ψ.sites)
nsite(ψ::zTMPS)     = length(ψ.sites_main)

siteindices(ψ::SignalMPS) = ψ.sites
siteindices(ψ::zTMPS)     = vcat(ψ.sites_main, ψ.sites_copy)
bondindices(ψ::SignalMPS) = ψ.bonds
bondindices(ψ::zTMPS)     = vcat(ψ.bonds_main, ψ.bonds_copy)

main_sites(ψ::zTMPS) = ψ.sites_main
copy_sites(ψ::zTMPS) = ψ.sites_copy
main_bonds(ψ::zTMPS) = ψ.bonds_main
copy_bonds(ψ::zTMPS) = ψ.bonds_copy

##################################### MPS UTILS #########################################

function as_signal_2n(ψ::zTMPS{I}) where {I<:Index}
    N = length(ψ.data)
    data  = Vector{ITensor{I}}(undef, 2N)
    sites = Vector{Index{I}}(undef, 2N)
    bonds = Vector{Index{I}}(undef, 2N - 1)

    for i in 1:N
        core = ψ.data[i]

        data[2i - 1]  = core.Amain
        data[2i]      = core.Acopy
        sites[2i - 1] = ψ.sites_main[i]
        sites[2i]     = ψ.sites_copy[i]

        # intra bond main(i)–copy(i)
        bonds[2i - 1] = ψ.bonds_copy[i]
        # inter bond copy(i)–main(i+1)
        if i < N
            bonds[2i] = ψ.bonds_main[i]
        end
    end

    return SignalMPS{I}(data, sites, bonds)
end

function writeback_signal_2n!(ψ::zTMPS{I}, ψ2n::SignalMPS{I}) where {I<:Index}
    N = length(ψ.data)
    length(ψ2n.sites) == 2N ||
        throw(ArgumentError("writeback_signal_2n!: 2n-sites mismatch."))

    for i in 1:N
        # Keep ψ.bonds_copy / ψ.bonds_main as the canonical bonds;
        # only overwrite the site tensors.
        ψ.data[i].Amain = ψ2n.data[2i - 1]
        ψ.data[i].Acopy = ψ2n.data[2i]
    end
    return ψ
end

end # module MPS
