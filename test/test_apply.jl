import QILaplace.ApplyMPO: apply

@testset "ApplyMPO.jl: Basic apply functionality" begin
    # Test 1: Apply identity MPO to a simple MPS
    @testset "Identity MPO application" begin
        n = 3
        # Create shared site indices
        sites = [Index(2, "site-$i") for i in 1:n]
        mpo_bonds = [Index(1, "mpo-bond-$i") for i in 1:(n - 1)]
        mps_bonds = [Index(2, "mps-bond-$i") for i in 1:(n - 1)]

        ψ = SignalMPS(sites, mps_bonds)

        # Create identity MPO with same sites
        W_data = Vector{ITensor}(undef, n)
        for i in 1:n
            s = sites[i]
            sp = prime(s)
            if i == 1
                W_data[i] = delta(sp, s) * ITensor([1.0], mpo_bonds[i])
            elseif i == n
                W_data[i] = delta(sp, s) * ITensor([1.0], mpo_bonds[i - 1])
            else
                W_data[i] =
                    delta(sp, s) *
                    ITensor([1.0], mpo_bonds[i - 1]) *
                    ITensor([1.0], mpo_bonds[i])
            end
        end
        W = SingleSiteMPO(W_data, sites, mpo_bonds)

        # Set MPS to random values
        for i in 1:n
            randn!(ψ.data[i])
        end

        # Apply identity MPO
        ψ_result = apply(W, ψ)

        # Check that result has correct structure
        @test length(ψ_result) == n
        @test length(ψ_result.sites) == n
        @test length(ψ_result.bonds) == n - 1

        # Check indices match
        @test ψ_result.sites == ψ.sites
    end

    # Test 2: Apply identity and compare dense vs apply result
    @testset "Identity MPO: dense vs apply" begin
        n = 4
        # Create shared site indices
        sites = [Index(2, "site-$i") for i in 1:n]
        mps_bonds = [Index(2, "mps-bond-$i") for i in 1:(n - 1)]
        mpo_bonds = [Index(1, "mpo-bond-$i") for i in 1:(n - 1)]

        ψ = SignalMPS(sites, mps_bonds)

        # Create identity MPO with same sites
        W_data = Vector{ITensor}(undef, n)
        for i in 1:n
            s = sites[i]
            sp = prime(s)
            if i == 1
                W_data[i] = delta(sp, s) * ITensor([1.0], mpo_bonds[i])
            elseif i == n
                W_data[i] = delta(sp, s) * ITensor([1.0], mpo_bonds[i - 1])
            else
                W_data[i] =
                    delta(sp, s) *
                    ITensor([1.0], mpo_bonds[i - 1]) *
                    ITensor([1.0], mpo_bonds[i])
            end
        end
        W = SingleSiteMPO(W_data, sites, mpo_bonds)

        # Initialize MPS with random values
        for i in 1:n
            randn!(ψ.data[i])
        end

        # Apply using the apply function
        ψ_apply = apply(W, ψ)

        # Apply using dense contraction
        result_apply_dense = to_dense_mps(ψ_apply)
        result_direct_dense = apply_dense(W, ψ)

        # Compare: they should be equal
        diff_tensor = array(result_apply_dense - result_direct_dense)
        @test LinearAlgebra.norm(diff_tensor) < 1e-12
    end
end

@testset "ApplyMPO.jl: Non-trivial MPO" begin
    @testset "Single-site gate MPO" begin
        n = 3

        # Create MPS and MPO with compatible indices
        sites = [Index(2, "site-$i") for i in 1:n]
        bonds = [Index(2, "bond-$i") for i in 1:(n - 1)]

        ψ = SignalMPS(sites, bonds)

        # Initialize MPS with random values
        for i in 1:n
            randn!(ψ.data[i])
        end

        # Create MPO with single-site Pauli-X gate on site 2
        W_data = Vector{ITensor}(undef, n)
        W_bonds = [Index(1, "W-bond-$i") for i in 1:(n - 1)]

        # Pauli X matrix
        X = [0.0 1.0; 1.0 0.0]

        for i in 1:n
            s = sites[i]
            sp = prime(s)

            if i == 1
                if n == 1
                    # Single site case
                    if i == 2  # Apply X gate
                        W_data[i] = ITensor(X, s', s)
                    else  # Identity
                        W_data[i] = delta(s', s)
                    end
                else
                    # First site with bond
                    if i == 2  # Apply X gate
                        W_data[i] = ITensor(X, s', s) * ITensor([1.0], W_bonds[i])
                    else  # Identity
                        W_data[i] = delta(s', s) * ITensor([1.0], W_bonds[i])
                    end
                end
            elseif i == n
                # Last site
                if i == 2  # Apply X gate
                    W_data[i] = ITensor(X, s', s) * ITensor([1.0], W_bonds[i - 1])
                else  # Identity
                    W_data[i] = delta(s', s) * ITensor([1.0], W_bonds[i - 1])
                end
            else
                # Bulk sites
                if i == 2  # Apply X gate
                    W_data[i] =
                        ITensor(X, s', s) *
                        ITensor([1.0], W_bonds[i - 1]) *
                        ITensor([1.0], W_bonds[i])
                else  # Identity
                    W_data[i] =
                        delta(s', s) *
                        ITensor([1.0], W_bonds[i - 1]) *
                        ITensor([1.0], W_bonds[i])
                end
            end
        end

        W = SingleSiteMPO(W_data, sites, W_bonds)

        # Apply using the apply function
        ψ_apply = apply(W, ψ)

        # Apply using dense contraction
        result_apply_dense = to_dense_mps(ψ_apply)
        result_direct_dense = apply_dense(W, ψ)

        # Compare
        @test LinearAlgebra.norm(result_apply_dense - result_direct_dense) < 1e-12
    end
end

@testset "ApplyMPO.jl: Random MPO application" begin
    @testset "Random MPO with small bond dimension" begin
        n = 4

        # Create compatible sites and bonds
        sites = [Index(2, "site-$i") for i in 1:n]
        mps_bonds = [Index(3, "mps-bond-$i") for i in 1:(n - 1)]
        mpo_bonds = [Index(2, "mpo-bond-$i") for i in 1:(n - 1)]

        # Create and initialize MPS
        ψ = SignalMPS(sites, mps_bonds)
        for i in 1:n
            randn!(ψ.data[i])
        end

        # Create and initialize MPO with random entries
        W_data = Vector{ITensor}(undef, n)
        for i in 1:n
            s = sites[i]
            sp = prime(s)

            if i == 1
                W_data[i] = random_itensor(sp, s, mpo_bonds[i])
            elseif i == n
                W_data[i] = random_itensor(mpo_bonds[i - 1], sp, s)
            else
                W_data[i] = random_itensor(mpo_bonds[i - 1], sp, s, mpo_bonds[i])
            end
        end

        W = SingleSiteMPO(W_data, sites, mpo_bonds)

        # Apply using the apply function
        ψ_apply = apply(W, ψ)

        # Apply using dense contraction
        result_apply_dense = to_dense_mps(ψ_apply)
        result_direct_dense = apply_dense(W, ψ)

        # Compare
        diff_norm = LinearAlgebra.norm(result_apply_dense - result_direct_dense)
        @test diff_norm < 1e-10
    end
end

@testset "ApplyMPO.jl: Error handling" begin
    @testset "Length mismatch" begin
        ψ = SignalMPS(3)
        W = SingleSiteMPO(4)

        @test_throws ArgumentError apply(W, ψ)
    end

    @testset "Site index mismatch" begin
        sites1 = [Index(2, "site-$i") for i in 1:3]
        sites2 = [Index(2, "site-$i") for i in 4:6]
        bonds = [Index(2, "bond-$i") for i in 1:2]

        ψ = SignalMPS(sites1, bonds)
        W_data = Vector{ITensor}(undef, 3)
        W_bonds = [Index(1, "W-bond-$i") for i in 1:2]

        for i in 1:3
            s = sites2[i]
            sp = prime(s)
            if i == 1
                W_data[i] = delta(sp, s) * ITensor([1.0], W_bonds[i])
            elseif i == 3
                W_data[i] = delta(sp, s) * ITensor([1.0], W_bonds[i - 1])
            else
                W_data[i] =
                    delta(sp, s) *
                    ITensor([1.0], W_bonds[i - 1]) *
                    ITensor([1.0], W_bonds[i])
            end
        end

        W = SingleSiteMPO(W_data, sites2, W_bonds)

        @test_throws ArgumentError apply(W, ψ)
    end

    @testset "MPO overlap mismatch" begin
        sites_a = [Index(2, "A-site-$i") for i in 1:2]
        bonds_a = [Index(2, "A-bond-1")]
        data_a = ITensor[
            random_itensor(prime(sites_a[1]), sites_a[1], bonds_a[1]),
            random_itensor(bonds_a[1], prime(sites_a[2]), sites_a[2]),
        ]
        W_a = SingleSiteMPO(data_a, sites_a, bonds_a)

        sites_b = [Index(2, "B-site-$i") for i in 1:2]
        bonds_b = [Index(2, "B-bond-1")]
        data_b = ITensor[
            random_itensor(prime(sites_b[1]), sites_b[1], bonds_b[1]),
            random_itensor(bonds_b[1], prime(sites_b[2]), sites_b[2]),
        ]
        W_b = SingleSiteMPO(data_b, sites_b, bonds_b)

        @test_throws ArgumentError apply(W_a, W_b)
    end
end

@testset "ApplyMPO.jl: PairedSiteMPO application" begin
    @testset "Identity PairedSiteMPO" begin
        n = 3
        W = PairedSiteMPO(n)
        ψ = zTMPS(W.sites_main, W.sites_copy)

        # Randomize ψ
        for core in ψ.data
            randn!(core.Amain)
            randn!(core.Acopy)
        end

        ψ_out = apply(W, ψ)

        # Compare using dense conversion
        ψ_2n = QILaplace.Mps._as_signal_2n(ψ)
        ψ_out_2n = QILaplace.Mps._as_signal_2n(ψ_out)

        v = to_dense_mps(ψ_2n)
        v_out = to_dense_mps(ψ_out_2n)

        @test LinearAlgebra.norm(v - v_out) < 1e-12
    end
end

@testset "ApplyMPO.jl: MPO composition" begin
    @testset "SingleSiteMPO composition via apply" begin
        n = 3
        sites = [Index(2, "site-$i") for i in 1:n]

        bonds1 = [Index(2, "b1-$i") for i in 1:(n - 1)]
        bonds2 = [Index(3, "b2-$i") for i in 1:(n - 1)]

        data1 = Vector{ITensor}(undef, n)
        data2 = Vector{ITensor}(undef, n)

        for i in 1:n
            s = sites[i]
            sp = prime(s)

            if i == 1
                data1[i] = random_itensor(sp, s, bonds1[i])
                data2[i] = random_itensor(sp, s, bonds2[i])
            elseif i == n
                data1[i] = random_itensor(bonds1[i - 1], sp, s)
                data2[i] = random_itensor(bonds2[i - 1], sp, s)
            else
                data1[i] = random_itensor(bonds1[i - 1], sp, s, bonds1[i])
                data2[i] = random_itensor(bonds2[i - 1], sp, s, bonds2[i])
            end
        end

        W1 = SingleSiteMPO(data1, sites, bonds1)
        W2 = SingleSiteMPO(data2, sites, bonds2)

        mps_bonds = [Index(2, "mps-bond-$i") for i in 1:(n - 1)]
        ψ = SignalMPS(sites, mps_bonds)
        for A in ψ.data
            randn!(A)
        end

        ψ_seq = apply(W2, apply(W1, ψ))
        W_comp = apply(W1, W2)
        ψ_comp = apply(W_comp, ψ)

        diff = LinearAlgebra.norm(array(to_dense_mps(ψ_seq) - to_dense_mps(ψ_comp)))
        @test diff < 1e-10

        # Test the contraction by explicitly multiplying dense tensors
        dense_direct = dense_compose_mpos(W1, W2)
        dense_comp = to_dense_mpo(W_comp)
        dense_direct = permute(dense_direct, inds(dense_comp))
        @test LinearAlgebra.norm(array(dense_direct - dense_comp)) < 1e-10
    end

    @testset "PairedSiteMPO composition via apply" begin
        n = 2

        sites_main = [Index(2, "sm-$i") for i in 1:n]
        sites_copy = [Index(2, "sc-$i") for i in 1:n]

        bonds_main1 = [Index(2, "bm1-$i") for i in 1:(n - 1)]
        bonds_copy1 = [Index(3, "bc1-$i") for i in 1:n]

        data1 = ITensor[
            random_itensor(sites_main[1], prime(sites_main[1]), bonds_copy1[1]),
            random_itensor(
                bonds_copy1[1], sites_copy[1], prime(sites_copy[1]), bonds_main1[1]
            ),
            random_itensor(
                bonds_main1[1], sites_main[2], prime(sites_main[2]), bonds_copy1[2]
            ),
            random_itensor(bonds_copy1[2], sites_copy[2], prime(sites_copy[2])),
        ]

        bonds_main2 = [Index(3, "bm2-$i") for i in 1:(n - 1)]
        bonds_copy2 = [Index(2, "bc2-$i") for i in 1:n]

        data2 = ITensor[
            random_itensor(sites_main[1], prime(sites_main[1]), bonds_copy2[1]),
            random_itensor(
                bonds_copy2[1], sites_copy[1], prime(sites_copy[1]), bonds_main2[1]
            ),
            random_itensor(
                bonds_main2[1], sites_main[2], prime(sites_main[2]), bonds_copy2[2]
            ),
            random_itensor(bonds_copy2[2], sites_copy[2], prime(sites_copy[2])),
        ]

        W1 = PairedSiteMPO(data1, sites_main, sites_copy, bonds_main1, bonds_copy1)
        W2 = PairedSiteMPO(data2, sites_main, sites_copy, bonds_main2, bonds_copy2)

        ψ = zTMPS(sites_main, sites_copy)
        for core in ψ.data
            randn!(core.Amain)
            randn!(core.Acopy)
        end

        ψ_seq = apply(W2, apply(W1, ψ))
        W_comp = apply(W1, W2)
        ψ_comp = apply(W_comp, ψ)

        dense_seq = to_dense_mps(QILaplace.Mps._as_signal_2n(ψ_seq))
        dense_comp = to_dense_mps(QILaplace.Mps._as_signal_2n(ψ_comp))

        diff = LinearAlgebra.norm(array(dense_seq - dense_comp))
        @test diff < 1e-10

        # Test the contraction by explicitly multiplying dense tensors
        single_W1 = QILaplace.ApplyMPO._as_single_site_mpo(W1)
        single_W2 = QILaplace.ApplyMPO._as_single_site_mpo(W2)
        single_W_comp = QILaplace.ApplyMPO._as_single_site_mpo(W_comp)
        dense_direct = dense_compose_mpos(single_W1, single_W2)
        dense_comp = to_dense_mpo(single_W_comp)
        dense_direct = permute(dense_direct, inds(dense_comp))
        @test LinearAlgebra.norm(array(dense_direct - dense_comp)) < 1e-9
    end

    @testset "Unequal length SingleSiteMPO composition" begin
        # Create W1 with 2 sites and W2 with 4 sites
        # W1 sites match W2.sites[2:3]
        all_sites = [Index(2, "site-$i") for i in 1:4]

        # W1: 2 sites (matches positions 2-3 of W2)
        sites1 = all_sites[2:3]
        bonds1 = [Index(2, "b1-1")]
        data1 = [
            random_itensor(sites1[1]', sites1[1], bonds1[1]),
            random_itensor(bonds1[1], sites1[2]', sites1[2]),
        ]
        W1 = SingleSiteMPO(data1, sites1, bonds1)

        # W2: 4 sites
        sites2 = all_sites
        bonds2 = [Index(3, "b2-$i") for i in 1:3]
        data2 = [
            random_itensor(sites2[1]', sites2[1], bonds2[1]),
            random_itensor(bonds2[1], sites2[2]', sites2[2], bonds2[2]),
            random_itensor(bonds2[2], sites2[3]', sites2[3], bonds2[3]),
            random_itensor(bonds2[3], sites2[4]', sites2[4]),
        ]
        W2 = SingleSiteMPO(data2, sites2, bonds2)

        # Compose: result should have 4 sites
        W_comp = apply(W1, W2)
        @test length(W_comp) == 4
        @test W_comp.sites == sites2

        # Basic structure check: result should have proper bond dimensions
        @test length(W_comp.bonds) == 3

        # Test the contraction by explicitly multiplying dense tensors
        W1_embedded = embed_mpo(W1, sites2)
        dense_direct = dense_compose_mpos(W1_embedded, W2)
        dense_comp = to_dense_mpo(W_comp)
        dense_direct = permute(dense_direct, inds(dense_comp))
        @test LinearAlgebra.norm(array(dense_direct - dense_comp)) < 1e-9
    end
end

@testset "ApplyMPO.jl: Operator interface" begin
    @testset "SingleSiteMPO * SignalMPS" begin
        n = 3
        sites = [Index(2, "op-site-$i") for i in 1:n]
        mps_bonds = [Index(2, "op-mps-bond-$i") for i in 1:(n - 1)]
        mpo_bonds = [Index(3, "op-mpo-bond-$i") for i in 1:(n - 1)]

        ψ = SignalMPS(sites, mps_bonds)
        for A in ψ.data
            randn!(A)
        end

        W_data = Vector{ITensor}(undef, n)
        for i in 1:n
            s = sites[i]
            sp = prime(s)
            if i == 1
                W_data[i] = random_itensor(sp, s, mpo_bonds[i])
            elseif i == n
                W_data[i] = random_itensor(mpo_bonds[i - 1], sp, s)
            else
                W_data[i] = random_itensor(mpo_bonds[i - 1], sp, s, mpo_bonds[i])
            end
        end
        W = SingleSiteMPO(W_data, sites, mpo_bonds)

        ψ_apply = apply(W, ψ)
        ψ_star = W * ψ

        diff = LinearAlgebra.norm(array(to_dense_mps(ψ_apply) - to_dense_mps(ψ_star)))
        @test diff < 1e-10
    end

    @testset "SingleSiteMPO * SingleSiteMPO" begin
        n = 3
        sites = [Index(2, "op2-site-$i") for i in 1:n]
        bonds1 = [Index(2, "op2-b1-$i") for i in 1:(n - 1)]
        bonds2 = [Index(3, "op2-b2-$i") for i in 1:(n - 1)]

        data1 = Vector{ITensor}(undef, n)
        data2 = Vector{ITensor}(undef, n)
        for i in 1:n
            s = sites[i]
            sp = prime(s)
            if i == 1
                data1[i] = random_itensor(sp, s, bonds1[i])
                data2[i] = random_itensor(sp, s, bonds2[i])
            elseif i == n
                data1[i] = random_itensor(bonds1[i - 1], sp, s)
                data2[i] = random_itensor(bonds2[i - 1], sp, s)
            else
                data1[i] = random_itensor(bonds1[i - 1], sp, s, bonds1[i])
                data2[i] = random_itensor(bonds2[i - 1], sp, s, bonds2[i])
            end
        end

        W1 = SingleSiteMPO(data1, sites, bonds1)
        W2 = SingleSiteMPO(data2, sites, bonds2)

        dense_apply = to_dense_mpo(apply(W1, W2))
        dense_star = to_dense_mpo(W1 * W2)
        dense_apply = permute(dense_apply, inds(dense_star))

        @test LinearAlgebra.norm(array(dense_apply - dense_star)) < 1e-10
    end

    @testset "Operator interface errors" begin
        # MPO * MPS length mismatch
        mps_sites = [Index(2, "err-site-$i") for i in 1:3]
        mps_bonds = [Index(2, "err-mps-bond-$i") for i in 1:2]
        ψ = SignalMPS(mps_sites, mps_bonds)

        mpo_sites = [Index(2, "err-mpo-site-$i") for i in 1:4]
        mpo_bonds = [Index(2, "err-mpo-bond-$i") for i in 1:3]
        mpo_data = Vector{ITensor}(undef, 4)
        for i in 1:4
            s = mpo_sites[i]
            sp = prime(s)
            if i == 1
                mpo_data[i] = random_itensor(sp, s, mpo_bonds[i])
            elseif i == 4
                mpo_data[i] = random_itensor(mpo_bonds[i - 1], sp, s)
            else
                mpo_data[i] = random_itensor(mpo_bonds[i - 1], sp, s, mpo_bonds[i])
            end
        end
        W_long = SingleSiteMPO(mpo_data, mpo_sites, mpo_bonds)

        @test_throws ArgumentError W_long * ψ

        # MPO * MPO with disjoint supports
        sites_left = [Index(2, "err-left-$i") for i in 1:2]
        bonds_left = [Index(2, "err-left-bond")]
        data_left = ITensor[
            random_itensor(prime(sites_left[1]), sites_left[1], bonds_left[1]),
            random_itensor(bonds_left[1], prime(sites_left[2]), sites_left[2]),
        ]
        W_left = SingleSiteMPO(data_left, sites_left, bonds_left)

        sites_right = [Index(2, "err-right-$i") for i in 1:2]
        bonds_right = [Index(2, "err-right-bond")]
        data_right = ITensor[
            random_itensor(prime(sites_right[1]), sites_right[1], bonds_right[1]),
            random_itensor(bonds_right[1], prime(sites_right[2]), sites_right[2]),
        ]
        W_right = SingleSiteMPO(data_right, sites_right, bonds_right)

        @test_throws ArgumentError W_left * W_right
    end
end
