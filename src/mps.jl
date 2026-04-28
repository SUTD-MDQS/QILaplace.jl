# src/Mps.jl
# This module contains the structures and functions related to Matrix Product States (MPS) that will be used in the main simulation code.
module Mps

using ITensors, Printf
import Base: getindex
import LinearAlgebra: norm

export SignalMPS,
    ZTMPS,
    nsite,
    siteindices,
    bondindices,
    update_site!,
    update_bond!,
    canonicalize!,
    compress!,
    coefficient,
    mps_to_vector

##################################### MPS TYPES #########################################

abstract type AbstractMPS end

"""
    PairCore(Amain::ITensor, Acopy::ITensor, c::Index)

A paired tensor core used inside [`ZTMPS`](@ref QILaplace.Mps.ZTMPS). Each `PairCore` holds two ITensors
(`Amain` and `Acopy`) that share exactly one common index `c` (the intra-site bond).
When a `ZTMPS` is printed, each site shows its `PairCore` components.

# Fields
- `Amain::ITensor` — tensor on the "main" register.
- `Acopy::ITensor` — tensor on the "copy" register.
- `c::Index`       — the shared intra-site bond index.
"""
mutable struct PairCore
    Amain::ITensor
    Acopy::ITensor
    c::Index
    function PairCore(Amain::ITensor, Acopy::ITensor, c::Index)
        # Rigorous sanity: Amain & Acopy must share exactly ONE common index and it must be c
        idxsA = collect(inds(Amain))
        idxsB = collect(inds(Acopy))
        shared = [I for I in idxsA if I in idxsB]
        length(shared) == 1 || throw(
            ArgumentError(
                "PairCore: Amain/Acopy must share exactly one index; found $(length(shared)) shared indices: $(shared)",
            ),
        )
        (shared[1] === c) ||
            throw(ArgumentError("PairCore: shared index $(shared[1]) ≠ provided c $(c)"))
        new(Amain, Acopy, c)
    end
end

"""
    SignalMPS{I<:Index} <: AbstractMPS

A Matrix Product State representing a 1D signal on `n` qubit sites.

# Fields
- `data::Vector{ITensor}` — the MPS tensors (one per site).
- `sites::Vector{I}`      — physical site indices.
- `bonds::Vector{I}`      — virtual bond indices (length `n-1`).
- `amplitude::Float64`    — the original Euclidean norm of the signal before normalisation.

Use [`signal_mps`](@ref QILaplace.SignalConverters.signal_mps) to construct a `SignalMPS` from a dense vector.
"""
mutable struct SignalMPS{I<:Index} <: AbstractMPS
    data::Vector{ITensor}
    sites::Vector{I}
    bonds::Vector{I}
    amplitude::Float64
    function SignalMPS(data::Vector{ITensor}, sites::Vector{I}, bonds::Vector{I}; amplitude::Float64=1.0) where {I}
        check_mps(data, sites, bonds)
        new{I}(data, sites, bonds, amplitude)
    end
end

"""
    ZTMPS{I<:Index} <: AbstractMPS

A paired-register MPS for non-unitary transforms (Damping Transform, z-Transform).
Each site consists of a [`PairCore`](@ref QILaplace.Mps.PairCore) holding a "main" and a "copy" tensor
connected by an intra-site bond.

# Fields
- `data::Vector{PairCore}` — paired tensor cores (one per site).
- `bonds_main::Vector{I}`  — inter-site bonds (length `n-1`).
- `bonds_copy::Vector{I}`  — intra-site bonds (length `n`).
- `sites_main::Vector{I}`  — physical site indices on the main register.
- `sites_copy::Vector{I}`  — physical site indices on the copy register.
- `amplitude::Float64`     — the original Euclidean norm of the signal.

Use [`signal_ztmps`](@ref QILaplace.SignalConverters.signal_ztmps) to construct a `ZTMPS` from a dense vector.
"""
mutable struct ZTMPS{I<:Index} <: AbstractMPS
    data::Vector{PairCore}
    bonds_main::Vector{I}
    bonds_copy::Vector{I}
    sites_main::Vector{I}
    sites_copy::Vector{I}
    amplitude::Float64

    function ZTMPS(
        cores::Vector{PairCore},
        bonds_main::Vector{I},
        bonds_copy::Vector{I},
        sites_main::Vector{I},
        sites_copy::Vector{I};
        amplitude::Float64=1.0,
    ) where {I}
        check_ztmps(cores, bonds_main, bonds_copy, sites_main, sites_copy)
        new{I}(cores, bonds_main, bonds_copy, sites_main, sites_copy, amplitude)
    end
end

##################################### MPS CONSTRUCTORS #########################################

function SignalMPS(sites::Vector{I}, bonds::Vector{I}) where {I}
    n = length(sites)
    n ≥ 1 || throw(ArgumentError("SignalMPS: Need at least 1 site to construct MPS"))
    length(bonds) == n - 1 || throw(
        ArgumentError(
            "SignalMPS: Number of bonds must be one less than number of sites. Got $(length(bonds)) bonds for $n sites.",
        ),
    )
    data = Vector{ITensor}(undef, n)
    if n == 1
        data[1] = ITensor(sites[1])
    else
        data[1] = ITensor(sites[1], bonds[1])
        for i in 2:(n-1)
            data[i] = ITensor(bonds[i-1], sites[i], bonds[i])
        end
        data[n] = ITensor(bonds[n-1], sites[n])
    end
    SignalMPS(data, sites, bonds)
end

SignalMPS(n::Int) = begin
    sites = [Index(2, @sprintf("site-%d", i)) for i in 1:n]
    bonds = [Index(1, @sprintf("bond-%d", i)) for i in 1:(n-1)]
    SignalMPS(sites, bonds)
end

function ZTMPS(sites_main::Vector{I}, sites_copy::Vector{I}) where {I}
    n = length(sites_main)
    n == length(sites_copy) ||
        throw(ArgumentError("sites_main and sites_copy must have same length"))

    bonds_main = Vector{I}(undef, max(n - 1, 0))
    for i in eachindex(bonds_main)
        bonds_main[i] = Index(1, @sprintf("bond-main-%d", i))::I
    end

    bonds_copy = Vector{I}(undef, n)
    for i in 1:n
        bonds_copy[i] = Index(1, @sprintf("bond-copy-%d", i))::I
    end

    data = Vector{PairCore}(undef, n)

    for i in 1:n
        c = bonds_copy[i] # intra-site bond main(i) <-> copy(i)

        bL = i == 1 ? nothing : bonds_main[i-1]
        bR = i == n ? nothing : bonds_main[i]

        Amain = bL === nothing ? ITensor(sites_main[i], c) : ITensor(bL, sites_main[i], c)
        Acopy = bR === nothing ? ITensor(c, sites_copy[i]) : ITensor(c, sites_copy[i], bR)

        data[i] = PairCore(Amain, Acopy, c)
    end

    return ZTMPS(data, bonds_main, bonds_copy, sites_main, sites_copy)
end

ZTMPS(n::Int) = begin
    sites_main = [Index(2, @sprintf("site-main-%d", i)) for i in 1:n]
    sites_copy = [Index(2, @sprintf("site-copy-%d", i)) for i in 1:n]
    ZTMPS(sites_main, sites_copy)
end

##################################### MPS VALIDATION #########################################

function check_mps(data::Vector{ITensor}, sites::Vector{I}, bonds::Vector{I}) where {I}
    n = length(data)
    length(sites) == n || throw(
        ArgumentError(
            "SignalMPS: Data length must equal number of sites. Got $(length(sites)) sites for $n tensors.",
        ),
    )
    length(bonds) == max(n - 1, 0) || throw(
        ArgumentError(
            "SignalMPS: Number of bonds must be one less than number of sites. Got $(length(bonds)) bonds for $n sites.",
        ),
    )

    if n == 1
        length(inds(data[1])) == 1 || throw(
            ArgumentError(
                "SignalMPS: Data tensor for single-site MPS must have exactly 1 index. Found $(length(inds(data[1])))",
            ),
        )
        (sites[1] in inds(data[1])) ||
            throw(ArgumentError("SignalMPS: Site index missing in single-site tensor"))
        return nothing
    end

    # Edge ranks = 2, bulk ranks = 3
    (length(inds(data[1])) == 2) || throw(
        ArgumentError(
            "SignalMPS: Edge tensors must have exactly 2 indices. Found $(length(inds(data[1]))) at site 1",
        ),
    )
    (length(inds(data[n])) == 2) || throw(
        ArgumentError(
            "SignalMPS: Edge tensors must have exactly 2 indices. Found $(length(inds(data[n]))) at site $n",
        ),
    )
    for i in 2:(n-1)
        length(inds(data[i])) == 3 || throw(
            ArgumentError(
                "SignalMPS: Bulk tensors must have exactly 3 indices. Found $(length(inds(data[i]))) at site $i",
            ),
        )
    end

    # Site presence and uniqueness
    for i in 1:n
        (sites[i] in inds(data[i])) ||
            throw(ArgumentError("SignalMPS: Site index missing in tensor at site $i"))
    end
    length(unique(sites)) == n ||
        throw(ArgumentError("SignalMPS: Site indices must be unique"))

    # Bond wiring (open boundary)
    for i in 1:(n-1)
        shared = commonind(inds(data[i]), inds(data[i+1]))
        (shared === bonds[i]) ||
            throw(ArgumentError("SignalMPS: Bond mismatch between core at $i and $(i+1)"))
    end
    length(unique(bonds)) == length(bonds) ||
        throw(ArgumentError("SignalMPS: Bond indices must be unique"))
    return nothing
end

check_mps(ψ::SignalMPS) = check_mps(ψ.data, ψ.sites, ψ.bonds)

# Iteration protocol for AbstractMPS
Base.iterate(ψ::SignalMPS) = iterate(ψ.data)
Base.iterate(ψ::SignalMPS, s) = iterate(ψ.data, s)
Base.iterate(ψ::ZTMPS) = iterate(ψ.data)
Base.iterate(ψ::ZTMPS, s) = iterate(ψ.data, s)

function check_ztmps(
    data::Vector{PairCore},
    bonds_main::Vector{I},
    bonds_copy::Vector{I},
    sites_main::Vector{I},
    sites_copy::Vector{I},
) where {I}
    N = length(data)

    length(sites_main) == N || throw(ArgumentError("ZTMPS: sites_main must have length n"))
    length(sites_copy) == N || throw(ArgumentError("ZTMPS: sites_copy must have length n"))
    length(bonds_copy) == N ||
        throw(ArgumentError("ZTMPS: bonds_copy (intra bonds) must have length N"))
    length(bonds_main) == max(N - 1, 0) ||
        throw(ArgumentError("ZTMPS: bonds_main (inter bonds) must have length N-1"))

    # Optional uniqueness checks (can relax if needed)
    length(unique(sites_main)) == N ||
        throw(ArgumentError("ZTMPS: sites_main must be unique"))
    length(unique(sites_copy)) == N ||
        throw(ArgumentError("ZTMPS: sites_copy must be unique"))
    length(unique(bonds_copy)) == N ||
        throw(ArgumentError("ZTMPS: bonds_copy (intra) must be unique"))
    if N > 1
        length(unique(bonds_main)) == length(bonds_main) ||
            throw(ArgumentError("ZTMPS: bonds_main (inter) must be unique"))
    end

    for i in 1:N
        core = data[i]
        isA = inds(core.Amain)
        isB = inds(core.Acopy)

        # site indices present
        sites_main[i] in isA ||
            throw(ArgumentError("ZTMPS: s_main[$i] missing in Amain[$i]"))
        sites_copy[i] in isB ||
            throw(ArgumentError("ZTMPS: s_copy[$i] missing in Acopy[$i]"))

        # intra bond consistency
        bonds_copy[i] == core.c ||
            throw(ArgumentError("ZTMPS: bonds_copy[$i] ≠ core.c at site $i"))
        core.c in isA ||
            throw(ArgumentError("ZTMPS: intra bond c[$i] missing in Amain[$i]"))
        core.c in isB ||
            throw(ArgumentError("ZTMPS: intra bond c[$i] missing in Acopy[$i]"))

        # left inter-site bond: copy(i-1) ↔ main(i)
        if i > 1
            bL = bonds_main[i-1]
            bL in isA || throw(
                ArgumentError("ZTMPS: left inter bond b_main[$(i-1)] missing in Amain[$i]"),
            )
            bL in inds(data[i-1].Acopy) || throw(
                ArgumentError("ZTMPS: left inter bond b_main[$(i-1)] not in Acopy[$(i-1)]"),
            )
        end

        # right inter-site bond: copy(i) ↔ main(i+1)
        if i < N
            bR = bonds_main[i]
            bR in isB || throw(
                ArgumentError("ZTMPS: right inter bond b_main[$i] missing in Acopy[$i]")
            )
        end
    end

    return nothing
end

function check_ztmps(ψ::ZTMPS)
    check_ztmps(ψ.data, ψ.bonds_main, ψ.bonds_copy, ψ.sites_main, ψ.sites_copy)
end

##################################### LIGHTWEIGHT API HELPERS #########################################

@deprecate nsite(ψ::AbstractMPS) length(ψ)
Base.length(ψ::SignalMPS) = length(ψ.data)
Base.length(ψ::ZTMPS) = length(ψ.data)

"""
    siteindices(ψ::SignalMPS)
    siteindices(ψ::ZTMPS)

Return a named tuple `(main=..., copy=...)` containing the physical site
indices of the MPS. For `SignalMPS`, `copy` is an empty vector.
"""
siteindices(ψ::SignalMPS) = (main=ψ.sites, copy=Vector{eltype(ψ.sites)}())
siteindices(ψ::ZTMPS) = (main=ψ.sites_main, copy=ψ.sites_copy)

"""
    bondindices(ψ::SignalMPS)
    bondindices(ψ::ZTMPS)

Return a named tuple `(main=..., copy=...)` containing the virtual bond
indices of the MPS. For `SignalMPS`, `copy` is an empty vector.
"""
bondindices(ψ::SignalMPS) = (main=ψ.bonds, copy=Vector{eltype(ψ.bonds)}())
bondindices(ψ::ZTMPS) = (main=ψ.bonds_main, copy=ψ.bonds_copy)

#################################### DISPLAY UTILS #########################################
function Base.show(io::IO, ψ::SignalMPS)
    N = length(ψ)
    println(io, "SignalMPS with $N sites:")
    for n in 1:N
        t = ψ.data[n]
        idxs = collect(inds(t))
        parts = String[]
        for I in idxs
            d = dim(I)
            tg = try
                s = string(tags(I))
                isempty(s) ? nothing : s
            catch
                nothing
            end
            push!(parts, tg === nothing ? "dim=$d" : "dim=$d, tags=$tg")
        end
        println(io, "  Site $n: ", join(parts, " | "))
    end
end

function Base.show(io::IO, ψ::ZTMPS)
    N = length(ψ)
    println(io, "ZTMPS with $N sites:")
    for n in 1:N
        core = ψ.data[n]
        println(io, "  Site $n:")
        # Amain
        tA = core.Amain
        idxsA = collect(inds(tA))
        partsA = String[]
        for I in idxsA
            d = dim(I)
            tg = try
                s = string(tags(I))
                isempty(s) ? nothing : s
            catch
                nothing
            end
            push!(partsA, tg === nothing ? "dim=$d" : "dim=$d, tags=$tg")
        end
        println(io, "    Amain: ", join(partsA, " | "))
        # Acopy
        tB = core.Acopy
        idxsB = collect(inds(tB))
        partsB = String[]
        for I in idxsB
            d = dim(I)
            tg = try
                s = string(tags(I))
                isempty(s) ? nothing : s
            catch
                nothing
            end
            push!(partsB, tg === nothing ? "dim=$d" : "dim=$d, tags=$tg")
        end
        println(io, "    Acopy: ", join(partsB, " | "))
    end
end

##################################### MPS UTILS #########################################
# Internal function for converting ZTMPS to SignalMPS
function _as_signal_2n(ψ::ZTMPS{I}) where {I}
    N = length(ψ.data)
    data = Vector{ITensor}(undef, 2N)
    sites = Vector{I}(undef, 2N)
    bonds = Vector{I}(undef, 2N - 1)

    for i in 1:N
        core = ψ.data[i]

        data[2i-1] = core.Amain
        data[2i] = core.Acopy
        sites[2i-1] = ψ.sites_main[i]
        sites[2i] = ψ.sites_copy[i]

        # intra bond main(i)–copy(i)
        bonds[2i-1] = ψ.bonds_copy[i]
        # inter bond copy(i)–main(i+1)
        if i < N
            bonds[2i] = ψ.bonds_main[i]
        end
    end

    return SignalMPS(data, sites, bonds; amplitude=ψ.amplitude)
end

# Internal function for converting SignalMPS to ZTMPS
function _writeback_signal_2n(ψ2n::SignalMPS{I}) where {I}
    N = length(ψ2n) ÷ 2
    length(ψ2n.sites) == 2N ||
        throw(ArgumentError("_writeback_signal_2n: 2n-sites mismatch."))
    length(ψ2n.bonds) == 2N - 1 ||
        throw(ArgumentError("_writeback_signal_2n: 2n-bonds mismatch."))
    zt_paircore_arr = Vector{PairCore}(undef, N)
    zt_sites_main = Vector{I}(undef, N)
    zt_sites_copy = Vector{I}(undef, N)
    zt_bonds_main = Vector{I}(undef, max(N - 1, 0))
    zt_bonds_copy = Vector{I}(undef, N)

    for i in 1:N
        zt_paircore_arr[i] = PairCore(ψ2n.data[2i-1], ψ2n.data[2i], ψ2n.bonds[2i-1])
        zt_sites_main[i] = ψ2n.sites[2i-1]
        zt_sites_copy[i] = ψ2n.sites[2i]
        zt_bonds_copy[i] = ψ2n.bonds[2i-1]
        if i < N
            zt_bonds_main[i] = ψ2n.bonds[2i]
        end
    end
    return ZTMPS(
        zt_paircore_arr, zt_bonds_main, zt_bonds_copy, zt_sites_main, zt_sites_copy;
        amplitude=ψ2n.amplitude,
    )
end

##################################### INDEX UPDATE FUNCTIONS #########################################

"""
    update_site!(ψ::SignalMPS, old_site_index::Index, new_site_index::Index)
    update_site!(ψ::ZTMPS, old_site_index::Index, new_site_index::Index)

Replace `old_site_index` with `new_site_index` in the MPS `ψ` (in-place).
Both indices must have the same dimension.
"""
function update_site!(ψ::SignalMPS{I}, old_site_index::I, new_site_index::I) where {I<:Index}
    # find in declared sites first
    site_idx = findfirst(x -> x == old_site_index, ψ.sites)
    site_idx === nothing && throw(
        ArgumentError(
            "update_site!: old site index $(old_site_index) not found in $(ψ.sites)"
        ),
    )

    # type consistency: site vector has element type E
    E = typeof(ψ.sites[site_idx])
    (new_site_index isa E) ||
        throw(ArgumentError("update_site!: new site index must be of type $E"))

    # dimension check
    dim(old_site_index) == dim(new_site_index) || throw(
        ArgumentError(
            "update_site!: Site index dimension mismatch at site $site_idx. Got $(dim(new_site_index)), found $(dim(old_site_index))",
        ),
    )

    # replace in core and site list
    replaceinds!(ψ.data[site_idx], old_site_index => new_site_index)
    ψ.sites[site_idx] = new_site_index
    return ψ
end

"""
    update_bond!(ψ::SignalMPS, old_bond_index::Index, new_bond_index::Index)
    update_bond!(ψ::ZTMPS, old_bond_index::Index, new_bond_index::Index)

Replace `old_bond_index` with `new_bond_index` in the MPS `ψ` (in-place).
Both indices must have the same dimension.
"""
function update_bond!(ψ::SignalMPS, old_bond_index::I, new_bond_index::I) where {I<:Index}
    bond_idx = findfirst(x -> x == old_bond_index, ψ.bonds)
    bond_idx === nothing && throw(
        ArgumentError(
            "update_bond!: old bond index $old_bond_index not found in $(ψ.bonds)"
        ),
    )

    B = typeof(ψ.bonds[bond_idx])
    (new_bond_index isa B) ||
        throw(ArgumentError("update_bond!: new bond index must be of type $B"))

    # dimension check
    dim(old_bond_index) == dim(new_bond_index) || throw(
        ArgumentError(
            "update_bond!: Bond index dimension mismatch at bond $bond_idx. Got $(dim(new_bond_index)), found $(dim(old_bond_index))",
        ),
    )

    # replace in neighbouring cores
    replaceinds!(ψ.data[bond_idx], old_bond_index => new_bond_index)
    replaceinds!(ψ.data[bond_idx+1], old_bond_index => new_bond_index)
    ψ.bonds[bond_idx] = new_bond_index
    return ψ
end

function update_site!(ψ::ZTMPS, old_site_index::I, new_site_index::I) where {I<:Index}
    # search main sites
    m = findfirst(x -> x == old_site_index, ψ.sites_main)
    if m !== nothing
        E = typeof(ψ.sites_main[m])
        (new_site_index isa E) ||
            throw(ArgumentError("update_site!: new main site index must be of type $E"))
        dim(old_site_index) == dim(new_site_index) || throw(
            ArgumentError("update_site!: Site index dimension mismatch at main site $m")
        )
        replaceinds!(ψ.data[m].Amain, old_site_index => new_site_index)
        ψ.sites_main[m] = new_site_index
        return ψ
    end
    # search copy sites
    c = findfirst(x -> x == old_site_index, ψ.sites_copy)
    if c !== nothing
        E = typeof(ψ.sites_copy[c])
        (new_site_index isa E) ||
            throw(ArgumentError("update_site!: new copy site index must be of type $E"))
        dim(old_site_index) == dim(new_site_index) || throw(
            ArgumentError("update_site!: Site index dimension mismatch at copy site $c")
        )
        replaceinds!(ψ.data[c].Acopy, old_site_index => new_site_index)
        ψ.sites_copy[c] = new_site_index
        return ψ
    end
    throw(ArgumentError("update_site!: old site index not found in ZTMPS"))
end

function update_bond!(ψ::ZTMPS, old_bond_index::I, new_bond_index::I) where {I<:Index}
    # search main bonds (inter-site bonds): bonds_main[k] connects copy(k) and main(k+1)
    k = findfirst(x -> x == old_bond_index, ψ.bonds_main)
    if k !== nothing
        B = typeof(ψ.bonds_main[k])
        (new_bond_index isa B) ||
            throw(ArgumentError("update_bond!: new main bond index must be of type $B"))
        dim(old_bond_index) == dim(new_bond_index) ||
            throw(ArgumentError("update_bond!: Main bond dimension mismatch at bond $k"))
        # replace in copy(k).Acopy and main(k+1).Amain
        replaceinds!(ψ.data[k].Acopy, old_bond_index => new_bond_index)
        replaceinds!(ψ.data[k+1].Amain, old_bond_index => new_bond_index)
        ψ.bonds_main[k] = new_bond_index
        return ψ
    end

    # search copy bonds (intra-site bonds): bonds_copy[j] connects Amain(j) and Acopy(j)
    j = findfirst(x -> x == old_bond_index, ψ.bonds_copy)
    if j !== nothing
        B = typeof(ψ.bonds_copy[j])
        (new_bond_index isa B) ||
            throw(ArgumentError("update_bond!: new copy bond index must be of type $B"))
        dim(old_bond_index) == dim(new_bond_index) ||
            throw(ArgumentError("update_bond!: Copy bond dimension mismatch at bond $j"))
        replaceinds!(ψ.data[j].Amain, old_bond_index => new_bond_index)
        replaceinds!(ψ.data[j].Acopy, old_bond_index => new_bond_index)
        ψ.bonds_copy[j] = new_bond_index
        ψ.data[j].c = new_bond_index
        return ψ
    end

    throw(ArgumentError("update_bond!: old bond index not found in ZTMPS"))
end

##################################### COEFFICIENT ACCESSORS #########################################

@inline function _site_projector(site::Index, raw::Integer)
    d = dim(site)
    val = Int(raw)
    (0 ≤ val < d) || throw(ArgumentError("coefficient: bit value $raw outside [0,$(d-1)]"))
    return setelt(dag(site) => (val + 1))
end

function _parse_config_string(spec::AbstractString)
    stripped = strip(spec)
    stripped = strip(stripped, ['[', ']', '(', ')', '{', '}'])
    isempty(stripped) && throw(ArgumentError("coefficient: configuration string is empty"))
    if occursin(r"[,\s]", stripped)
        tokens = split(stripped, r"[,\s]+"; keepempty=false)
        isempty(tokens) && throw(
            ArgumentError("coefficient: configuration string did not contain any entries"),
        )
        return parse.(Int, tokens)
    else
        all(c -> c in ('0', '1'), stripped) ||
            throw(ArgumentError("coefficient: bit strings may contain only '0' or '1'"))
        return [c == '1' ? 1 : 0 for c in stripped]
    end
end

function _bits_from_integer(value::Integer, n::Int)
    value ≥ 0 ||
        throw(ArgumentError("coefficient: integer configuration must be non-negative"))
    bits = Vector{Int}(undef, n)
    tmp = value
    for i in n:-1:1
        bits[i] = Int(tmp & 0x1)
        tmp >>= 1
    end
    tmp == 0 ||
        throw(ArgumentError("coefficient: integer $value requires more than $n bits"))
    return bits
end

"""
    coefficient(ψ::SignalMPS, config)
    coefficient(ψ::ZTMPS, config)

Return the amplitude ⟨config|ψ⟩ for the given zero-based bit configuration.

`config` can be:
- A `Vector` or `Tuple` of integers, each in `[0, dim(site)-1]`.
- A bit string like `"1010"` or `"[1,0,1,0]"`.
- A non-negative integer interpreted as an `n`-bit big-endian pattern.

You may also use direct indexing: `ψ[0, 1, 0, 1]` is equivalent to
`coefficient(ψ, (0, 1, 0, 1))`.

# Examples
```julia
ψ = signal_mps([1.0, 0.0, 0.0, 0.0])
coefficient(ψ, [0, 0])   # amplitude of |00⟩
coefficient(ψ, "00")      # same, via bit string
ψ[0, 0]                   # same, via indexing
```
"""
function coefficient(ψ::SignalMPS, config::AbstractVector{<:Integer})
    N = length(ψ.data)
    length(config) == N ||
        throw(ArgumentError("coefficient: expected $N entries, got $(length(config))"))
    amp = ψ.data[1] * _site_projector(ψ.sites[1], config[1])
    @inbounds for i in 2:N
        amp = (amp * ψ.data[i]) * _site_projector(ψ.sites[i], config[i])
    end
    return ψ.amplitude * scalar(amp)
end

coefficient(ψ::SignalMPS, config::Tuple{Vararg{Integer}}) = coefficient(ψ, collect(config))
coefficient(ψ::SignalMPS, config::Vararg{Integer}) = coefficient(ψ, collect(config))
coefficient(ψ::SignalMPS, spec::AbstractString) = coefficient(ψ, _parse_config_string(spec))
function coefficient(ψ::SignalMPS, value::Integer)
    coefficient(ψ, _bits_from_integer(value, length(ψ.data)))
end

function coefficient(ψ::ZTMPS, config)
    ψ_2n = _as_signal_2n(ψ)
    return coefficient(ψ_2n, config)
end

getindex(ψ::SignalMPS, config::Vararg{Integer}) = coefficient(ψ, collect(config))
getindex(ψ::ZTMPS, config::Vararg{Integer}) = coefficient(ψ, collect(config))

##################################### MPS TO VECTOR #########################################

"""
    mps_to_vector(ψ::SignalMPS; reverse::Bool=false)

Extract the dense state vector from a `SignalMPS`, scaled by the stored `amplitude`.

# Arguments
- `reverse::Bool=false`: If `false` (default), return the vector in MSB-first
  (normal/natural) bit ordering — matching the original signal ordering from
  `signal_mps`. If `true`, return the vector in the raw column-major
  (bit-reversed) ordering of the tensor network, which corresponds to the
  output ordering after a QFT or other transform.

# Examples
```julia
ψ = signal_mps([1.0, 2.0, 3.0, 4.0])
mps_to_vector(ψ)                  # [1.0, 2.0, 3.0, 4.0] (original signal)
mps_to_vector(ψ; reverse=true)    # bit-reversed ordering
```
"""
function mps_to_vector(ψ::SignalMPS; reverse::Bool=false)
    n = length(ψ.sites)
    N = 2^n

    T = ψ.data[1]
    for i in 2:n
        T *= ψ.data[i]
    end

    site_order = reverse ? ψ.sites : Base.reverse(ψ.sites)
    arr = Array(T, site_order...)
    return reshape(arr, N) * ψ.amplitude
end

"""
    mps_to_vector(ψ::ZTMPS; reverse::Bool=false)

Extract the dense state vector from a `ZTMPS` by first converting to its
`2N`-site `SignalMPS` representation, then extracting as a flat vector.
Scaled by the stored `amplitude`.

# Arguments
- `reverse::Bool=false`: Controls bit ordering (see `mps_to_vector(::SignalMPS)`).
"""
function mps_to_vector(ψ::ZTMPS; reverse::Bool=false)
    ψ2n = _as_signal_2n(ψ)
    return mps_to_vector(ψ2n; reverse=reverse)
end

##################################### NORM FUNCTIONS #########################################

"""
    norm(ψ::SignalMPS)
    norm(ψ::ZTMPS)

Compute the Euclidean norm √⟨ψ|ψ⟩ of the MPS state by full tensor contraction.
This extends `LinearAlgebra.norm`.
"""
function norm(ψ::SignalMPS)
    # Create conjugate with primed bonds to prevent bond contraction
    ψ_conj = [dag(replaceinds(t, (b => prime(b) for b in ψ.bonds)...)) for t in ψ.data]

    # Site indices automatically match (same MPS), so they contract
    E = ITensor(1)
    @inbounds for i in 1:length(ψ.data)
        E *= ψ.data[i]
        E *= ψ_conj[i]
    end
    return sqrt(abs(scalar(E)))
end

# norm(ψ::ZTMPS) — delegates to the SignalMPS method via _as_signal_2n.
function norm(ψ::ZTMPS)
    ψ = _as_signal_2n(ψ)
    return norm(ψ)
end

##################################### MPS CANONICALIZATION #########################################

"""
    canonicalize!(ψ, direction; center=nothing, cutoff=1e-12, maxdim=typemax(Int))

Bring the MPS `ψ` into canonical form (in-place) via QR/LQ sweeps.

# Arguments
- `direction::String`: `"->"` (left-canonical) or `"<-"` (right-canonical).
- `center::Int`: orthogonality center site (defaults to `N` for `"->"`, `1` for `"<-"`).
- `cutoff`, `maxdim`: truncation parameters.

Works for both `SignalMPS` and `ZTMPS`.
"""
function canonicalize!(
    ψ::SignalMPS,
    direction::Symbol;
    center::Union{Nothing,Int}=nothing,
    cutoff::Float64=1e-12,
    maxdim::Int=typemax(Int),
)
    direction ∈ (:right, :left) || throw(ArgumentError("Direction must be :right or :left"))

    N = length(ψ.data)

    if direction == :right
        c = something(center, N)
        1 ≤ c ≤ N || throw(DomainError(c, "Center out of range [1,$N]"))

        for n in 1:(c-1)
            left_inds = n == 1 ? (ψ.sites[n],) : (ψ.bonds[n-1], ψ.sites[n])
            U, R = factorize(
                ψ.data[n], left_inds...; ortho="left", cutoff=cutoff, maxdim=maxdim
            )
            newlink = commonind(U, R)
            newlink !== nothing || throw(ErrorException("No common link at site $n"))

            canon = Index(dim(newlink); tags=tags(ψ.bonds[n]))
            replaceinds!(U, newlink => canon)
            replaceinds!(R, newlink => canon)

            ψ.data[n] = U
            ψ.data[n+1] = R * ψ.data[n+1]
            ψ.bonds[n] = canon
        end
    else  # direction == :left
        c = something(center, 1)
        1 ≤ c ≤ N || throw(DomainError(c, "Center out of range [1,$N]"))

        for n in N:-1:(c+1)
            left_inds = (ψ.bonds[n-1],)
            L, V = factorize(
                ψ.data[n], left_inds...; ortho="right", cutoff=cutoff, maxdim=maxdim
            )
            newlink = commonind(L, V)
            newlink !== nothing || throw(ErrorException("No common link at site $n"))

            canon = Index(dim(newlink); tags=tags(ψ.bonds[n-1]))
            replaceinds!(L, newlink => canon)
            replaceinds!(V, newlink => canon)

            ψ.data[n-1] = ψ.data[n-1] * L
            ψ.data[n] = V
            ψ.bonds[n-1] = canon
        end
    end

    check_mps(ψ)
    return ψ
end

# Convenience wrapper with automatic center selection
function canonicalize!(ψ::SignalMPS, direction::Symbol, cutoff::Float64, maxdim::Int)
    canonicalize!(ψ, direction; center=nothing, cutoff=cutoff, maxdim=maxdim)
end

"""
    canonicalize!(ψ::ZTMPS, direction::Symbol; center::Union{Nothing,Int}=nothing,
                  cutoff::Float64=1e-12, maxdim::Int=typemax(Int))

Bring ZTMPS into canonical form by sweeping through PairCore structures.

# Arguments
- `ψ`: ZTMPS to canonicalize in-place
- `direction`: `:right` (left-to-right) or `:left` (right-to-left)
- `center`: orthogonality center pair index (defaults to n for `:right`, 1 for `:left`)
- `cutoff`: truncation threshold for singular values
- `maxdim`: maximum bond dimension

# Details
For each pair, canonicalizes both Amain and Acopy tensors while maintaining
the PairCore structure. Sweeps through pairs in the specified direction.
"""
function canonicalize!(
    ψ::ZTMPS,
    direction::Symbol;
    center::Union{Nothing,Int}=nothing,
    cutoff::Float64=1e-12,
    maxdim::Int=typemax(Int),
)
    direction ∈ (:right, :left) || throw(ArgumentError("Direction must be :right or :left"))

    # Convert to 2N SignalMPS, canonicalize there, and write back into ψ (mutate in-place)
    ψ_2n = _as_signal_2n(ψ)
    canonicalize!(ψ_2n, direction; center=center, cutoff=cutoff, maxdim=maxdim)
    ψ_back = _writeback_signal_2n(ψ_2n)

    # Mutate ψ in-place to reflect ψ_back (preserve same object identity)
    N = length(ψ.data)
    for i in 1:N
        ψ.data[i].Amain = ψ_back.data[i].Amain
        ψ.data[i].Acopy = ψ_back.data[i].Acopy
        ψ.data[i].c = ψ_back.data[i].c
        ψ.sites_main[i] = ψ_back.sites_main[i]
        ψ.sites_copy[i] = ψ_back.sites_copy[i]
        ψ.bonds_copy[i] = ψ_back.bonds_copy[i]
        if i < N
            ψ.bonds_main[i] = ψ_back.bonds_main[i]
        end
    end

    check_ztmps(ψ.data, ψ.bonds_main, ψ.bonds_copy, ψ.sites_main, ψ.sites_copy)
    return ψ
end

# Convenience wrapper with automatic center selection
function canonicalize!(ψ::ZTMPS, direction::Symbol, cutoff::Float64, maxdim::Int)
    canonicalize!(ψ, direction; center=nothing, cutoff=cutoff, maxdim=maxdim)
end

##################################### MPS COMPRESSION #########################################

"""
    compress!(ψ; maxdim=typemax(Int), tol=1e-12, sweeps=1)

Truncate bond dimensions of the MPS `ψ` via alternating SVD sweeps (in-place).
After compression, the MPS is re-normalised to unit norm.

Works for both `SignalMPS` and `ZTMPS`.
"""
function compress!(
    ψ::SignalMPS{I}; maxdim::Int=typemax(Int), tol::Float64=1e-12, sweeps::Int=1
) where {I}
    # SVD-based two-site compression swept left-right and right-left
    N = length(ψ.data)
    N ≥ 2 || throw(DomainError("SignalMPS must have at least 2 sites."))

    cutoff = tol^2 / ((N - 1) * sweeps)

    # Make right-orthogonal to start
    canonicalize!(ψ, :left)

    for _ in 1:sweeps
        # Left -> Right
        for j in 1:(N-1)
            left_inds = j == 1 ? (ψ.sites[j],) : (ψ.sites[j], ψ.bonds[j-1])
            U, S, V = svd(
                ψ.data[j] * ψ.data[j+1], left_inds...; cutoff=cutoff, maxdim=maxdim
            )
            svd_link = commonind(U, S)
            new_bond = Index(dim(svd_link); tags=tags(ψ.bonds[j]))

            U_left, U_right = U, S * V
            replaceinds!(U_left, svd_link => new_bond)
            replaceinds!(U_right, svd_link => new_bond)
            ψ.bonds[j] = new_bond

            ψ.data[j] = U_left
            ψ.data[j+1] = U_right
        end
        # Right -> Left
        for j in (N-1):-1:1
            left_inds = j == 1 ? (ψ.sites[j],) : (ψ.bonds[j-1], ψ.sites[j])
            U, S, V = svd(
                ψ.data[j] * ψ.data[j+1], left_inds...; cutoff=cutoff, maxdim=maxdim
            )
            svd_link = commonind(S, V)
            new_bond = Index(dim(svd_link); tags=tags(ψ.bonds[j]))

            U_left, U_right = U * S, V
            replaceinds!(U_left, svd_link => new_bond)
            replaceinds!(U_right, svd_link => new_bond)
            ψ.bonds[j] = new_bond

            ψ.data[j] = U_left
            ψ.data[j+1] = U_right
        end
    end

    # Re-canonicalize and normalize
    canonicalize!(ψ, :left)
    check_mps(ψ)

    # Normalize: absorb norm into amplitude, keep tensor data unit-normed
    nrm = norm(ψ)
    if nrm != 0
        ψ.amplitude *= nrm
        ψ.data[1] *= 1.0 / nrm
    end
    return ψ
end

function compress!(
    ψ::ZTMPS{I}; maxdim::Int=typemax(Int), tol::Float64=1e-12, sweeps::Int=1
) where {I}
    # Convert to 2N SignalMPS, compress, then write back in-place
    ψ2n = _as_signal_2n(ψ)
    compress!(ψ2n; maxdim=maxdim, tol=tol, sweeps=sweeps)
    ψ_back = _writeback_signal_2n(ψ2n)

    # Mutate input ZTMPS to reflect compressed tensors
    N = length(ψ.data)
    for i in 1:N
        ψ.data[i].Amain = ψ_back.data[i].Amain
        ψ.data[i].Acopy = ψ_back.data[i].Acopy
        ψ.data[i].c = ψ_back.data[i].c
        ψ.sites_main[i] = ψ_back.sites_main[i]
        ψ.sites_copy[i] = ψ_back.sites_copy[i]
        ψ.bonds_copy[i] = ψ_back.bonds_copy[i]
        if i < N
            ψ.bonds_main[i] = ψ_back.bonds_main[i]
        end
    end

    check_ztmps(ψ.data, ψ.bonds_main, ψ.bonds_copy, ψ.sites_main, ψ.sites_copy)
    return ψ
end

end # module MPS
