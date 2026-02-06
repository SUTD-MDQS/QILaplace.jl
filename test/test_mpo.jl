import QILaplace.Mpo: SingleSiteMPO, PairedSiteMPO, check_singlesitempo, check_pairedsitempo
import QILaplace.Mps: PairCore

@testset "mpo.jl: SingleSiteMPO struct validation" begin
    # valid SingleSiteMPO for N = 1, 2, 3, 4 (must produce the given mpos without errors)
    for N in 1:4
        W = SingleSiteMPO(N)
        @test length(W.data) == N
        @test length(W.sites) == N
        @test length(W.bonds) == max(N-1, 0)
    end

    # SingleSiteMpo must have data, sites, bonds fields
    W = SingleSiteMPO(2)
    @test hasfield(typeof(W), :data)
    @test hasfield(typeof(W), :sites)
    @test hasfield(typeof(W), :bonds)

    # SingleSiteMPO(N) must return delta function MPO at each site
    W = SingleSiteMPO(3)
    for i in 1:3
        t = W.data[i]
        s = W.sites[i]
        sp = s'
        # Contract out dummy bond indices (dimension 1) to get delta(s, s')
        t_contracted = copy(t)
        for idx in inds(t)
            if idx != s && idx != sp
                t_contracted *= ITensor(1.0, idx)
            end
        end
        # Check diagonal elements are 1
        for d in 1:dim(s)
            val = t_contracted[s => d, sp => d]
            @test abs(val - 1.0) < 1e-10
        end
    end

    # length of sites must equal length of data. Must error otherwise
    s = [Index(2, "s1"), Index(2, "s2")]
    b = [Index(1, "b1")]
    d = [random_itensor(s[1], s[1]', b[1])]  # only 1 tensor for 2 sites
    @test_throws ArgumentError SingleSiteMPO(d, s, b)

    # length of bonds must equal length of data - 1. Must error otherwise
    s = [Index(2, "s1"), Index(2, "s2")]
    b = [Index(1, "b1"), Index(1, "b2")]  # 2 bonds for 2 sites
    d = [random_itensor(s[1], s[1]', b[1]), random_itensor(s[2], s[2]', b[1])]
    @test_throws ArgumentError SingleSiteMPO(d, s, b)

    # For n=1, must have empty bonds and single site index with rank 2
    s = [Index(2, "s1")]
    d = [random_itensor(s[1], s[1]')]
    W = SingleSiteMPO(d, s, eltype(s)[])
    @test length(W.bonds) == 0
    @test length(inds(W.data[1])) == 2

    # Error if n=1 has wrong rank
    s = [Index(2, "s1")]
    b_fake = Index(1, "b_fake")
    d = [random_itensor(s[1], s[1]', b_fake)]  # rank 3 instead of 2
    @test_throws ArgumentError SingleSiteMPO(d, s, eltype(s)[])

    # Edge rank must be 3 and bulk rank must be 4 for n >= 2
    N = 4
    s = [Index(2, @sprintf("s%d", i)) for i in 1:N]
    b = [Index(1, @sprintf("b%d", i)) for i in 1:(N - 1)]
    d = [
        random_itensor(s[1], s[1]', b[1]),
        random_itensor(b[1], s[2], s[2]', b[2]),
        random_itensor(b[2], s[3], s[3]', b[3]),
        random_itensor(b[3], s[4], s[4]'),
    ]
    W = SingleSiteMPO(d, s, b)
    @test length(inds(W.data[1])) == 3  # edge
    @test length(inds(W.data[2])) == 4  # bulk
    @test length(inds(W.data[3])) == 4  # bulk
    @test length(inds(W.data[4])) == 3  # edge

    # Error if edge rank is wrong
    bad_d = [
        random_itensor(s[1], s[1]'),  # rank 2 instead of 3
        random_itensor(b[1], s[2], s[2]', b[2]),
        random_itensor(b[2], s[3], s[3]'),
    ]
    @test_throws ArgumentError SingleSiteMPO(bad_d, s[1:3], b[1:2])

    # In each site, the site and its prime must be present
    for i in 1:N
        @test s[i] in inds(W.data[i])
        @test s[i]' in inds(W.data[i])
    end

    # Error if site index is missing
    bad_d = [
        random_itensor(b[1]),  # missing site index
        random_itensor(b[1], s[2], s[2]', b[2]),
        random_itensor(b[2], s[3], s[3]'),
    ]
    @test_throws ArgumentError SingleSiteMPO(bad_d, s[1:3], b[1:2])

    # Site indices must be unique
    dup_s = [s[1], s[1], s[3], s[4]]
    @test_throws ArgumentError SingleSiteMPO(d, dup_s, b)

    # Bonds must be unique indices
    dup_b = [b[1], b[1], b[1]]
    @test_throws ArgumentError SingleSiteMPO(d, s, dup_b)

    # bonds must exist in adjacent data cores
    wrong_bond = Index(1, "wrong")
    bad_b = [b[1], wrong_bond, b[3]]
    @test_throws ArgumentError SingleSiteMPO(d, s, bad_b)

    # check_singlesitempo must work on a SingleSiteMPO argument
    W = SingleSiteMPO(3)
    @test check_singlesitempo(W) === nothing
    @test check_singlesitempo(W.data, W.sites, W.bonds) === nothing
end

@testset "mpo.jl: PairedSiteMPO struct validation" begin
    # valid PairedSiteMPO for N = 1, 2, 3, 4 (must produce the given mpos without errors)
    for N in 1:4
        W = PairedSiteMPO(N)
        @test length(W.data) == 2N
        @test length(W.sites_main) == N
        @test length(W.sites_copy) == N
        @test length(W.bonds_main) == max(N-1, 0)
        @test length(W.bonds_copy) == N
    end

    # PairedSiteMpo must have data, sites_main, sites_copy, bonds_main, bonds_copy fields
    W = PairedSiteMPO(2)
    @test hasfield(typeof(W), :data)
    @test hasfield(typeof(W), :sites_main)
    @test hasfield(typeof(W), :sites_copy)
    @test hasfield(typeof(W), :bonds_main)
    @test hasfield(typeof(W), :bonds_copy)

    # PairedSiteMPO(N) must return delta function MPO at each site
    W = PairedSiteMPO(3)
    for i in 1:3
        # Check main sites - contract out dummy bonds and verify diagonal
        s_main = W.sites_main[i]
        t_main = copy(W.data[2i - 1])
        for idx in inds(W.data[2i - 1])
            if idx != s_main && idx != s_main'
                t_main *= ITensor(1.0, idx)
            end
        end
        for d in 1:dim(s_main)
            val = t_main[s_main => d, s_main' => d]
            @test abs(val - 1.0) < 1e-10
        end

        # Check copy sites - contract out dummy bonds and verify diagonal
        s_copy = W.sites_copy[i]
        t_copy = copy(W.data[2i])
        for idx in inds(W.data[2i])
            if idx != s_copy && idx != s_copy'
                t_copy *= ITensor(1.0, idx)
            end
        end
        for d in 1:dim(s_copy)
            val = t_copy[s_copy => d, s_copy' => d]
            @test abs(val - 1.0) < 1e-10
        end
    end

    # length of sites_main + sites_copy must equal length of data respectively
    n = 2
    sm = [Index(2, "sm1"), Index(2, "sm2")]
    sc = [Index(2, "sc1"), Index(2, "sc2")]
    bm = [Index(1, "bm1")]
    bc = [Index(1, "bc1"), Index(1, "bc2")]
    d = [random_itensor(sm[1], sm[1]', bc[1])]  # only 1 tensor instead of 4
    @test_throws ArgumentError PairedSiteMPO(d, sm, sc, bm, bc)

    # length of bonds_main must be N-1 and bonds_copy must be N
    n = 2
    sm = [Index(2, "sm1"), Index(2, "sm2")]
    sc = [Index(2, "sc1"), Index(2, "sc2")]
    bm = [Index(1, "bm1"), Index(1, "bm2")]  # should be 1, not 2
    bc = [Index(1, "bc1"), Index(1, "bc2")]
    d = [
        random_itensor(sm[1], sm[1]', bc[1]),
        random_itensor(sc[1], sc[1]', bc[1], bm[1]),
        random_itensor(sm[2], sm[2]', bm[1], bc[2]),
        random_itensor(sc[2], sc[2]', bc[2]),
    ]
    @test_throws ArgumentError PairedSiteMPO(d, sm, sc, bm, bc)

    # For n=1, both data must have rank 3, common bond must be the bond_copy
    n = 1
    sm = [Index(2, "sm1")]
    sc = [Index(2, "sc1")]
    bm = eltype(sm)[]
    bc = [Index(1, "bc1")]
    d = [random_itensor(sm[1], sm[1]', bc[1]), random_itensor(sc[1], sc[1]', bc[1])]
    W = PairedSiteMPO(d, sm, sc, bm, bc)
    @test length(inds(W.data[1])) == 3
    @test length(inds(W.data[2])) == 3
    @test commonind(W.data[1], W.data[2]) === bc[1]
    @test length(W.bonds_main) == 0
    @test length(W.bonds_copy) == 1

    # Error if n=1 has wrong bond
    bad_bc = [Index(1, "bad")]
    d_bad = [random_itensor(sm[1], sm[1]', bc[1]), random_itensor(sc[1], sc[1]', bc[1])]
    @test_throws ArgumentError PairedSiteMPO(d_bad, sm, sc, bm, bad_bc)

    # Edge rank must be 3 and bulk rank must be 4 for n >= 2
    n = 3
    sm = [Index(2, @sprintf("sm%d", i)) for i in 1:n]
    sc = [Index(2, @sprintf("sc%d", i)) for i in 1:n]
    bm = [Index(1, @sprintf("bm%d", i)) for i in 1:(n - 1)]
    bc = [Index(1, @sprintf("bc%d", i)) for i in 1:n]
    d = [
        random_itensor(sm[1], sm[1]', bc[1]),
        random_itensor(sc[1], sc[1]', bc[1], bm[1]),
        random_itensor(sm[2], sm[2]', bm[1], bc[2]),
        random_itensor(sc[2], sc[2]', bc[2], bm[2]),
        random_itensor(sm[3], sm[3]', bm[2], bc[3]),
        random_itensor(sc[3], sc[3]', bc[3]),
    ]
    W = PairedSiteMPO(d, sm, sc, bm, bc)
    @test length(inds(W.data[1])) == 3   # edge (main1)
    @test length(inds(W.data[2])) == 4   # bulk (copy1)
    @test length(inds(W.data[3])) == 4   # bulk (main2)
    @test length(inds(W.data[4])) == 4   # bulk (copy2)
    @test length(inds(W.data[5])) == 4   # bulk (main3)
    @test length(inds(W.data[6])) == 3   # edge (copy3)

    # In each site (main and copy), the site and its prime must be present
    for i in 1:n
        @test sm[i] in inds(W.data[2i - 1])
        @test sm[i]' in inds(W.data[2i - 1])
        @test sc[i] in inds(W.data[2i])
        @test sc[i]' in inds(W.data[2i])
    end

    # Error if site index is missing
    bad_d = [
        random_itensor(bc[1]),  # missing site index
        random_itensor(sc[1], sc[1]', bc[1]),
    ]
    @test_throws ArgumentError PairedSiteMPO(bad_d, sm[1:1], sc[1:1], bm, bc[1:1])

    # Site indices (main and copy) must be unique
    dup_sm = [sm[1], sm[1], sm[3]]
    @test_throws ArgumentError PairedSiteMPO(d, dup_sm, sc, bm, bc)

    dup_sc = [sc[1], sc[1], sc[3]]
    @test_throws ArgumentError PairedSiteMPO(d, sm, dup_sc, bm, bc)

    # Bonds (main and copy) must be unique indices
    dup_bm = [bm[1], bm[1]]
    @test_throws ArgumentError PairedSiteMPO(d, sm, sc, dup_bm, bc)

    dup_bc = [bc[1], bc[1], bc[3]]
    @test_throws ArgumentError PairedSiteMPO(d, sm, sc, bm, dup_bc)

    # bonds must exist in adjacent data cores (main and copy)
    wrong_bond = Index(1, "wrong")
    bad_bm = [wrong_bond, bm[2]]
    @test_throws ArgumentError PairedSiteMPO(d, sm, sc, bad_bm, bc)

    bad_bc = [bc[1], wrong_bond, bc[3]]
    @test_throws ArgumentError PairedSiteMPO(d, sm, sc, bm, bad_bc)

    # check_pairedsitempo must work on a PairedSiteMPO argument
    W = PairedSiteMPO(3)
    @test check_pairedsitempo(W) === nothing
    @test check_pairedsitempo(
        W.data, W.sites_main, W.sites_copy, W.bonds_main, W.bonds_copy
    ) === nothing
end

@testset "mpo.jl: Base.show methods" begin
    # Must print things as expected for SingleSiteMPO
    W_single = SingleSiteMPO(2)
    io = IOBuffer()
    show(io, W_single)
    output = String(take!(io))
    @test occursin("SingleSiteMPO with 2 sites", output)
    @test occursin("Site 1:", output)
    @test occursin("Site 2:", output)

    # Must print things as expected for PairedSiteMPO
    W_paired = PairedSiteMPO(2)
    io = IOBuffer()
    show(io, W_paired)
    output = String(take!(io))
    @test occursin("PairedSiteMPO with 2 sites", output)
    @test occursin("Site 1 (main):", output)
    @test occursin("Site 1 (copy):", output)
    @test occursin("Site 2 (main):", output)
    @test occursin("Site 2 (copy):", output)
end

@testset "mpo.jl: SingleSiteMPO site and bond index update" begin
    N = 3
    W = SingleSiteMPO(N)

    # update_site! must replace the old site index with the new site index in the sites array
    old_site = W.sites[2]
    new_site = Index(dim(old_site), "site-updated-2")
    QILaplace.Mpo.update_site!(W, old_site, new_site)
    @test W.sites[2] == new_site
    @test new_site in inds(W.data[2])
    @test new_site' in inds(W.data[2])
    # Verify structure still valid
    @test check_singlesitempo(W) === nothing

    # update_bond! must replace the old bond index with the new bond index in the bonds array
    old_bond = W.bonds[1]
    new_bond = Index(dim(old_bond), "bond-updated-1")
    QILaplace.Mpo.update_bond!(W, old_bond, new_bond)
    @test W.bonds[1] == new_bond
    @test new_bond in inds(W.data[1])
    @test new_bond in inds(W.data[2])
    # Verify structure still valid
    @test check_singlesitempo(W) === nothing

    # must throw error if a site to be updated has dimension mismatch with the new site
    bad_site = Index(dim(W.sites[1]) + 1, "bad-dim-site")
    @test_throws ArgumentError QILaplace.Mpo.update_site!(W, W.sites[1], bad_site)

    # must throw error if dimension mismatch for bond
    bad_bond = Index(dim(W.bonds[1]) + 1, "bad-dim-bond")
    @test_throws ArgumentError QILaplace.Mpo.update_bond!(W, W.bonds[1], bad_bond)

    # must throw error if a site is not found in the sites array
    fake_site = Index(2, "fake-site")
    new_fake = Index(2, "new-fake")
    @test_throws ArgumentError QILaplace.Mpo.update_site!(W, fake_site, new_fake)

    # must throw error if a bond is not found in the bonds array
    fake_bond = Index(1, "fake-bond")
    new_fake_bond = Index(1, "new-fake-bond")
    @test_throws ArgumentError QILaplace.Mpo.update_bond!(W, fake_bond, new_fake_bond)
end

@testset "mpo.jl: PairedSiteMPO site and bond index update" begin
    n = 3
    W = PairedSiteMPO(n)

    # update_site! must replace the old site index in sites_main
    old_site_main = W.sites_main[2]
    new_site_main = Index(dim(old_site_main), "smain-updated-2")
    QILaplace.Mpo.update_site!(W, old_site_main, new_site_main)
    @test W.sites_main[2] == new_site_main
    @test new_site_main in inds(W.data[2 * 2 - 1])
    @test new_site_main' in inds(W.data[2 * 2 - 1])
    # Verify structure still valid
    @test check_pairedsitempo(W) === nothing

    # update_site! must replace the old site index in sites_copy
    old_site_copy = W.sites_copy[1]
    new_site_copy = Index(dim(old_site_copy), "scopy-updated-1")
    QILaplace.Mpo.update_site!(W, old_site_copy, new_site_copy)
    @test W.sites_copy[1] == new_site_copy
    @test new_site_copy in inds(W.data[2 * 1])
    @test new_site_copy' in inds(W.data[2 * 1])
    # Verify structure still valid
    @test check_pairedsitempo(W) === nothing

    # update_bond! must replace the old bond index in bonds_main
    old_bond_main = W.bonds_main[1]
    new_bond_main = Index(dim(old_bond_main), "bmain-updated-1")
    QILaplace.Mpo.update_bond!(W, old_bond_main, new_bond_main)
    @test W.bonds_main[1] == new_bond_main
    @test new_bond_main in inds(W.data[2 * 1])      # copy1
    @test new_bond_main in inds(W.data[2 * 1 + 1])  # main2
    # Verify structure still valid
    @test check_pairedsitempo(W) === nothing

    # update_bond! must replace the old bond index in bonds_copy
    old_bond_copy = W.bonds_copy[2]
    new_bond_copy = Index(dim(old_bond_copy), "bcopy-updated-2")
    QILaplace.Mpo.update_bond!(W, old_bond_copy, new_bond_copy)
    @test W.bonds_copy[2] == new_bond_copy
    @test new_bond_copy in inds(W.data[2 * 2 - 1])  # main2
    @test new_bond_copy in inds(W.data[2 * 2])      # copy2
    # Verify structure still valid
    @test check_pairedsitempo(W) === nothing

    # must throw error if a site to be updated has dimension mismatch
    bad_site_main = Index(dim(W.sites_main[1]) + 1, "bad-dim-main")
    @test_throws ArgumentError QILaplace.Mpo.update_site!(W, W.sites_main[1], bad_site_main)

    bad_site_copy = Index(dim(W.sites_copy[1]) + 1, "bad-dim-copy")
    @test_throws ArgumentError QILaplace.Mpo.update_site!(W, W.sites_copy[1], bad_site_copy)

    # must throw error if dimension mismatch for bonds
    bad_bond_main = Index(dim(W.bonds_main[1]) + 1, "bad-dim-bmain")
    @test_throws ArgumentError QILaplace.Mpo.update_bond!(W, W.bonds_main[1], bad_bond_main)

    bad_bond_copy = Index(dim(W.bonds_copy[1]) + 1, "bad-dim-bcopy")
    @test_throws ArgumentError QILaplace.Mpo.update_bond!(W, W.bonds_copy[1], bad_bond_copy)

    # must throw error if a site is not found in sites_main/sites_copy
    fake_site = Index(2, "fake-site")
    new_fake = Index(2, "new-fake")
    @test_throws ArgumentError QILaplace.Mpo.update_site!(W, fake_site, new_fake)

    # must throw error if a bond is not found in bonds_main/bonds_copy
    fake_bond = Index(1, "fake-bond")
    new_fake_bond = Index(1, "new-fake-bond")
    @test_throws ArgumentError QILaplace.Mpo.update_bond!(W, fake_bond, new_fake_bond)
end
