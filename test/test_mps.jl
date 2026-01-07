using Test, QILaplace, ITensors, Printf
using LinearAlgebra

import QILaplace.Mps: _as_signal_2n, _writeback_signal_2n, canonicalize!, compress!, norm

# PairCore mutable struct check
@testset "mps.jl: PairCore struct validation" begin
    # Valid PairCore with exactly one common index
    c = Index(2, "common")
    sA = Index(2, "siteA")
    sB = Index(2, "siteB")
    Amain = random_itensor(sA, c)
    Acopy = random_itensor(c, sB)
    
    pc = QILaplace.Mps.PairCore(Amain, Acopy, c)
    @test pc.Amain === Amain
    @test pc.Acopy === Acopy
    @test pc.c === c
    
    # Error when more than 1 common bond
    c2 = Index(2, "common2")
    Amain_bad = random_itensor(sA, c, c2)
    Acopy_bad = random_itensor(c, c2, sB)
    @test_throws ArgumentError QILaplace.Mps.PairCore(Amain_bad, Acopy_bad, c)
    
    # Error when no common bond
    Amain_none = random_itensor(sA)
    Acopy_none = random_itensor(sB)
    c_fake = Index(2, "fake")
    @test_throws ArgumentError QILaplace.Mps.PairCore(Amain_none, Acopy_none, c_fake)
end

# SignalMPS mutable struct check
@testset "mps.jl: SignalMPS struct validation" begin
    # Valid SignalMPS for N = 1, 2, 3, 4
    for N in 1:4
        ψ = QILaplace.Mps.SignalMPS(N)
        @test length(ψ.data) == N
        @test length(ψ.sites) == N
        @test length(ψ.bonds) == max(N-1, 0)
    end
    
    # Length mismatch: data ≠ sites
    s = [Index(2, "s1"), Index(2, "s2")]
    b = [Index(2, "b1")]
    d = [random_itensor(s[1], b[1])]  # only 1 tensor for 2 sites
    @test_throws ArgumentError QILaplace.Mps.SignalMPS(d, s, b)
    
    # Edge ranks = 2, bulk = 3
    N = 3
    sites = [Index(2, @sprintf("s%d", i)) for i in 1:N]
    bonds = [Index(2, @sprintf("b%d", i)) for i in 1:N-1]
    data = [random_itensor(sites[1], bonds[1]),
            random_itensor(bonds[1], sites[2], bonds[2]),
            random_itensor(bonds[2], sites[3])]
    ψ = QILaplace.Mps.SignalMPS(data, sites, bonds)
    @test length(inds(ψ.data[1])) == 2
    @test length(inds(ψ.data[2])) == 3
    @test length(inds(ψ.data[3])) == 2
    
    # Site index missing
    bad_data = [random_itensor(bonds[1]),  # missing site index
                random_itensor(bonds[1], sites[2], bonds[2]),
                random_itensor(bonds[2], sites[3])]
    @test_throws ArgumentError QILaplace.Mps.SignalMPS(bad_data, sites, bonds)
    
    # Non-unique sites
    dup_sites = [sites[1], sites[1], sites[3]]
    @test_throws ArgumentError QILaplace.Mps.SignalMPS(data, dup_sites, bonds)
    
    # Bond mismatch
    wrong_bond = Index(2, "wrong")
    bad_bonds = [bonds[1], wrong_bond]
    @test_throws ArgumentError QILaplace.Mps.SignalMPS(data, sites, bad_bonds)
end

# zTMPS mutable struct check
@testset "mps.jl: zTMPS struct validation" begin
    # Valid zTMPS for n = 1, 2, 3, 4
    for n in 1:4
        ψ = QILaplace.Mps.zTMPS(n)
        @test length(ψ.data) == n
        @test length(ψ.sites_main) == n
        @test length(ψ.sites_copy) == n
        @test length(ψ.bonds_main) == max(n-1, 0)
        @test length(ψ.bonds_copy) == n
        
        for i in 1:n
            @test ψ.data[i] isa QILaplace.Mps.PairCore
            @test ψ.sites_main[i] in inds(ψ.data[i].Amain)
            @test ψ.sites_copy[i] in inds(ψ.data[i].Acopy)
            @test ψ.bonds_copy[i] in inds(ψ.data[i].Amain)
            @test ψ.bonds_copy[i] in inds(ψ.data[i].Acopy)
            @test ψ.bonds_copy[i] == ψ.data[i].c
        end
        
        if n > 1
            for i in 1:n-1
                @test ψ.bonds_main[i] in inds(ψ.data[i].Acopy)
                @test ψ.bonds_main[i] in inds(ψ.data[i+1].Amain)
            end
        end
        
        # Site uniqueness check
        @test length(unique([ψ.sites_main; ψ.sites_copy])) == 2*nsite(ψ)
        
        # Bond uniqueness check
        @test length(unique([ψ.bonds_main; ψ.bonds_copy])) == length(ψ.bonds_main) + length(ψ.bonds_copy)
    end
end

# Base.show method tests
@testset "mps.jl: Base.show methods" begin
    # SignalMPS show output
    ψs = QILaplace.Mps.SignalMPS(2)
    io = IOBuffer()
    show(io, ψs)
    output = String(take!(io))
    @test occursin("SignalMPS with 2 sites", output)
    @test occursin("Site 1:", output)
    @test occursin("Site 2:", output)
    
    # zTMPS show output
    ψz = QILaplace.Mps.zTMPS(2)
    io = IOBuffer()
    show(io, ψz)
    output = String(take!(io))
    @test occursin("zTMPS with 2 sites", output)
    @test occursin("Amain:", output)
    @test occursin("Acopy:", output)
end

# Round-trip internal conversion between SignalMPS and zTMPS
@testset "mps.jl: as_signal_2n / writeback_signal_2n round-trip" begin
    n = 3
    ψ = QILaplace.Mps.zTMPS(n)
    for i in 1:n
        ψ.data[i].Amain .= random_itensor(inds(ψ.data[i].Amain)...) 
        ψ.data[i].Acopy .= random_itensor(inds(ψ.data[i].Acopy)...)
    end

    ψ2n = _as_signal_2n(ψ)
    @test length(ψ2n.data) == 2n

    ψ_back = _writeback_signal_2n(ψ2n)
    for i in 1:n
        @test inds(ψ.data[i].Amain) == inds(ψ_back.data[i].Amain)
        @test inds(ψ.data[i].Acopy) == inds(ψ_back.data[i].Acopy)
        @test ψ_back.data[i].c == ψ.data[i].c
    end
end

# Canonicalization via conversion
@testset "mps.jl: canonicalize via conversion" begin
    n = 4
    ψ = QILaplace.Mps.zTMPS(n)
    for i in 1:n
        ψ.data[i].Amain .= random_itensor(inds(ψ.data[i].Amain)...) 
        ψ.data[i].Acopy .= random_itensor(inds(ψ.data[i].Acopy)...) 
    end

    ψ2n = _as_signal_2n(ψ)
    before = norm(ψ)

    canonicalize!(ψ2n, "->")
    canonicalize!(ψ2n, "<-")

    ψ_back = _writeback_signal_2n(ψ2n)
    after = norm(ψ_back)
    @test isapprox(before, after; rtol=1e-10)
    QILaplace.Mps.check_ztmps(ψ_back.data, ψ_back.bonds_main, ψ_back.bonds_copy, ψ_back.sites_main, ψ_back.sites_copy)
end

# update_site! and update_bonds! index updates for SignalMPS
@testset "mps.jl: SignalMPS site and bond index update" begin
    # SignalMPS updates should yield a valid MPS
    N = 4
    sites = [Index(2, @sprintf("s%d", i)) for i in 1:N]
    bonds = [Index(2, @sprintf("b%d", i)) for i in 1:N-1]
    data = [random_itensor(sites[1], bonds[1]),
            random_itensor(bonds[1], sites[2], bonds[2]),
            random_itensor(bonds[2], sites[3], bonds[3]),
            random_itensor(bonds[3], sites[4])]
    ψs = QILaplace.Mps.SignalMPS(data, sites, bonds)

    # update site index 2
    old_site = ψs.sites[2]
    new_site = Index(dim(old_site), "site-updated-2")
    QILaplace.Mps.update_site!(ψs, old_site, new_site)
    # After update we should still pass MPS structural validation
    QILaplace.Mps.check_mps(ψs)
    @test ψs.sites[2] == new_site

    # update bond 2
    old_bond = ψs.bonds[2]
    new_bond = Index(dim(old_bond), "bond-updated-2")
    QILaplace.Mps.update_bond!(ψs, old_bond, new_bond)
    QILaplace.Mps.check_mps(ψs)
    @test ψs.bonds[2] == new_bond
end

# update_site! and update_bonds! index updates for zTMPS
@testset "mps.jl: zTMPS site and bond index update" begin
    n = 3
    ψz = QILaplace.Mps.zTMPS(n)
    for i in 1:n
        ψz.data[i].Amain .= random_itensor(inds(ψz.data[i].Amain)...)
        ψz.data[i].Acopy .= random_itensor(inds(ψz.data[i].Acopy)...)
    end

    # update a main site
    old_sm = ψz.sites_main[2]
    new_sm = Index(dim(old_sm), "smain-upd-2")
    QILaplace.Mps.update_site!(ψz, old_sm, new_sm)
    QILaplace.Mps.check_ztmps(ψz.data, ψz.bonds_main, ψz.bonds_copy, ψz.sites_main, ψz.sites_copy)
    @test ψz.sites_main[2] == new_sm

    # update a copy site
    old_sc = ψz.sites_copy[1]
    new_sc = Index(dim(old_sc), "scopy-upd-1")
    QILaplace.Mps.update_site!(ψz, old_sc, new_sc)
    QILaplace.Mps.check_ztmps(ψz.data, ψz.bonds_main, ψz.bonds_copy, ψz.sites_main, ψz.sites_copy)
    @test ψz.sites_copy[1] == new_sc

    # update a main bond
    old_bm = ψz.bonds_main[1]
    new_bm = Index(dim(old_bm), "bmain-upd-1")
    QILaplace.Mps.update_bond!(ψz, old_bm, new_bm)
    QILaplace.Mps.check_ztmps(ψz.data, ψz.bonds_main, ψz.bonds_copy, ψz.sites_main, ψz.sites_copy)

    # update a copy bond
    old_bc = ψz.bonds_copy[2]
    new_bc = Index(dim(old_bc), "bcopy-upd-2")
    QILaplace.Mps.update_bond!(ψz, old_bc, new_bc)
    QILaplace.Mps.check_ztmps(ψz.data, ψz.bonds_main, ψz.bonds_copy, ψz.sites_main, ψz.sites_copy)
    @test ψz.data[2].c == new_bc

    # Dimension mismatch on bond update should throw
    bad_old = ψz.bonds_copy[1]
    bad_new = Index(dim(bad_old) + 1, "bad-dim-bond")
    @test_throws ArgumentError QILaplace.Mps.update_bond!(ψz, bad_old, bad_new)

    # Site update with mismatched dimension should throw
    bad_old_site = ψz.sites_main[1]
    bad_new_site = Index(dim(bad_old_site) + 1, "bad-dim-site")
    @test_throws ArgumentError QILaplace.Mps.update_site!(ψz, bad_old_site, bad_new_site)
end

# norm function verification with analytical norm for both SignalMPS and zTMPS
@testset "mps.jl: norm function verification" begin
    # SignalMPS norm - contract to array and compare with LinearAlgebra.norm
    N = 3
    sites = [Index(2, @sprintf("s%d", i)) for i in 1:N]
    bonds = [Index(2, @sprintf("b%d", i)) for i in 1:N-1]
    data = [random_itensor(sites[1], bonds[1]),
            random_itensor(bonds[1], sites[2], bonds[2]),
            random_itensor(bonds[2], sites[3])]
    ψs = QILaplace.Mps.SignalMPS(data, sites, bonds)
    
    # Contract MPS to dense array
    ψs_contracted = ψs.data[1]
    for i in 2:N
        ψs_contracted *= ψs.data[i]
    end
    ψs_array = Array(ψs_contracted, sites...)
    analytical_norm = LinearAlgebra.norm(ψs_array)
    
    n1 = norm(ψs)
    @test isapprox(n1, analytical_norm, rtol=1e-10)
    @test n1 > 0
    @test isfinite(n1)
    
    # zTMPS norm - contract to array and compare
    ψz = QILaplace.Mps.zTMPS(3)
    for i in 1:3
        ψz.data[i].Amain .= random_itensor(inds(ψz.data[i].Amain)...)
        ψz.data[i].Acopy .= random_itensor(inds(ψz.data[i].Acopy)...)
    end
    
    # Contract zTMPS via conversion to SignalMPS
    ψz_2n = _as_signal_2n(ψz)
    ψz_contracted = ψz_2n.data[1]
    for i in 2:length(ψz_2n.data)
        ψz_contracted *= ψz_2n.data[i]
    end
    ψz_array = Array(ψz_contracted, ψz_2n.sites...)
    analytical_norm_z = LinearAlgebra.norm(ψz_array)
    
    n2 = norm(ψz)
    @test isapprox(n2, analytical_norm_z, rtol=1e-10)
    @test n2 > 0
    @test isfinite(n2)
    
    # Norm scaling - verify against analytical value
    ψz.data[1].Amain *= 2.0
    n3 = norm(ψz)
    @test isapprox(n3, 2.0 * n2, rtol=1e-10)
    
    # Verify scaling with analytical norm
    ψz_2n_scaled = _as_signal_2n(ψz)
    ψz_contracted_scaled = ψz_2n_scaled.data[1]
    for i in 2:length(ψz_2n_scaled.data)
        ψz_contracted_scaled *= ψz_2n_scaled.data[i]
    end
    ψz_array_scaled = Array(ψz_contracted_scaled, ψz_2n_scaled.sites...)
    analytical_norm_z_scaled = LinearAlgebra.norm(ψz_array_scaled)
    @test isapprox(n3, analytical_norm_z_scaled, rtol=1e-10)
end

# Compression tests
@testset "mps.jl: compression and norm preservation" begin
    N = 4
    sites = [Index(2, @sprintf("s%d", i)) for i in 1:N]
    bonds = [Index(4, @sprintf("b%d", i)) for i in 1:N-1]
    data = [random_itensor(sites[1], bonds[1]),
            random_itensor(bonds[1], sites[2], bonds[2]),
            random_itensor(bonds[2], sites[3], bonds[3]),
            random_itensor(bonds[3], sites[4])]
    ψ = QILaplace.Mps.SignalMPS(data, sites, bonds)

    QILaplace.Mps.compress!(ψ; maxdim=2, tol=1e-8, sweeps=2)
    n_after = norm(ψ)
    @test isapprox(n_after, 1.0; rtol=1e-8)
    for b in ψ.bonds
        @test dim(b) ≤ 2
    end

    # zTMPS compress via conversion
    n = 3
    sm = [Index(2, @sprintf("sm%d", i)) for i in 1:n]
    sc = [Index(2, @sprintf("sc%d", i)) for i in 1:n]
    ψz = QILaplace.Mps.zTMPS(sm, sc)
    for i in 1:n
        ψz.data[i].Amain .= random_itensor(inds(ψz.data[i].Amain)...)
        ψz.data[i].Acopy .= random_itensor(inds(ψz.data[i].Acopy)...)
    end
    QILaplace.Mps.compress!(ψz; maxdim=2, tol=1e-8, sweeps=2)
    @test isapprox(norm(ψz), 1.0; rtol=1e-8)
    for b in ψz.bonds_main; @test dim(b) ≤ 2; end
    for b in ψz.bonds_copy; @test dim(b) ≤ 2; end
end

# Additional tests: alternation invariants and writeback validation
@testset "mps.jl: alternation invariants and writeback validation" begin
    n = 3
    ψ = QILaplace.Mps.zTMPS(n)
    for i in 1:n
        ψ.data[i].Amain .= random_itensor(inds(ψ.data[i].Amain)...)
        ψ.data[i].Acopy .= random_itensor(inds(ψ.data[i].Acopy)...)
    end

    ψ2n = _as_signal_2n(ψ)
    # Check alternating site/bond ordering
    for i in 1:n
        @test ψ2n.data[2i - 1] === ψ.data[i].Amain
        @test ψ2n.data[2i]     === ψ.data[i].Acopy
        @test ψ2n.sites[2i - 1] == ψ.sites_main[i]
        @test ψ2n.sites[2i]     == ψ.sites_copy[i]
        @test ψ2n.bonds[2i - 1] == ψ.bonds_copy[i]
        if i < n
            @test ψ2n.bonds[2i] == ψ.bonds_main[i]
        end
    end

    # Corrupt the ordering (swap two bonds) and ensure writeback rejects it
    bad = deepcopy(ψ2n)
    bad.bonds[1], bad.bonds[2] = bad.bonds[2], bad.bonds[1]
    @test_throws ArgumentError _writeback_signal_2n(bad)

    # Corrupt the sites (swap a main/copy pair) and ensure writeback rejects it
    bad2 = deepcopy(ψ2n)
    bad2.sites[1], bad2.sites[2] = bad2.sites[2], bad2.sites[1]
    @test_throws ArgumentError _writeback_signal_2n(bad2)
end
