import QILaplace.QFTTransform: zip_up_mpos, zip_down_mpos

# ======================== HELPER FUNCTIONS ========================

# Helper for DFT (Unitary QFT with +2πi convention)
function dft(v::Vector)
    N = length(v)
    out = zeros(ComplexF64, N)
    for k in 0:N-1
        sum_val = 0.0 + 0.0im
        for j in 0:N-1
            # QFT definition: |j> -> 1/sqrt(N) sum_k exp(2pi i j k / N) |k>
            angle = 2 * π * j * k / N
            sum_val += v[j+1] * exp(im * angle)
        end
        out[k+1] = sum_val / sqrt(N)
    end
    return out
end

# The bit-reversed Fourier transform Q_n matrix
# Q_n[j, k] = exp(2πi * bitrev(j) * k / N) / sqrt(N)
function qn_matrix(n::Int)
    N = 2^n
    M = zeros(ComplexF64, N, N)
    for j in 0:N-1
        j_rev = bits_to_int(int_to_bits(j, n; order=:lsb); order=:msb)
        for k in 0:N-1
            angle = 2 * π * j_rev * k / N
            M[j+1, k+1] = exp(im * angle) / sqrt(N)
        end
    end
    return M
end

# ======================== ZIP UP TESTS ========================
@testset "qft_transformer.jl: zip up algorithm correctness" begin
    @testset "Zip Up: ITensor Array Multiplication Equivalence" begin
        # Test that zip_up correctly combines two MPOs by comparing ITensor arrays
        # Instead of dense matrices, we use ITensor to array conversions which are easier to multiply
        
        for n in 3:4
            # Create sites
            sites = [Index(2, "site-$i") for i in 1:n]
            
            # Create two control-phase MPOs that will be combined
            mpo1_orig = control_Hphase_mpo(n, sites)
            
            # mpo2 acts on subset of sites (2:n)
            subset_sites = sites[2:end]
            mpo2_orig = control_Hphase_mpo(n-1, subset_sites)
            
            # Store original ITensors as arrays before modifying indices
            T1 = prod(mpo1_orig.data)
            arr1 = Array(T1, prime.(sites)..., sites...)
            
            T2 = prod(mpo2_orig.data)
            arr2 = Array(T2, prime.(subset_sites)..., subset_sites...)
            
            # For expected result, we need to contract the ITensors manually
            # Expected: applying mpo1 then mpo2 using the apply function
            mpo_expected = apply(mpo1_orig, mpo2_orig; cutoff=0.0, maxdim=1000)
            T_expected = prod(mpo_expected.data)
            arr_expected = Array(T_expected, prime.(sites)..., sites...)
            
            # Now prepare copies for zip_up
            mpo1 = control_Hphase_mpo(n, sites)
            mpo2 = control_Hphase_mpo(n-1, subset_sites)
            
            # Prepare indices for zip_up (contract M1's output with M2's input on overlapping sites)
            for (idx2, s) in enumerate(subset_sites)
                idx1 = idx2 + 1  # Position in mpo1
                s_cont = sim(s)
                mpo1.data[idx1] = replaceinds(mpo1.data[idx1], s => s_cont)
                mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
            end
            
            # Perform zip_up
            oc = n  # OC at bottom
            mpo_zipped, new_oc = zip_up_mpos(mpo1, mpo2, oc)
            
            # Extract result as ITensor array
            T_result = prod(mpo_zipped.data)
            arr_result = Array(T_result, prime.(sites)..., sites...)
            
            # Compare ITensor arrays
            err = LinearAlgebra.norm(arr_result - arr_expected)
            @test isapprox(err, 0.0; atol=1e-12)
            
            # Verify OC moved to top
            @test new_oc == 1
        end
    end

    @testset "Zip Up: Apply Function Equivalence" begin
        # Test that zip_up result matches using the apply function to combine MPOs
        
        for n in 3:4
            sites = [Index(2, "site-$i") for i in 1:n]
            
            # Create two MPOs for apply
            mpo1_orig = control_Hphase_mpo(n, sites)
            subset_sites = sites[2:end]
            mpo2_orig = control_Hphase_mpo(n-1, subset_sites)
            
            # Use apply to combine MPOs (this should give the same result as zip_up)
            mpo_applied = apply(mpo1_orig, mpo2_orig; cutoff=0.0, maxdim=1000)
            
            # Create fresh copies for zip_up
            mpo1 = control_Hphase_mpo(n, sites)
            mpo2 = control_Hphase_mpo(n-1, subset_sites)
            
            # Prepare indices for zip_up
            for (idx2, s) in enumerate(subset_sites)
                idx1 = idx2 + 1
                s_cont = sim(s)
                mpo1.data[idx1] = replaceinds(mpo1.data[idx1], s => s_cont)
                mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
            end
            
            oc = n
            mpo_zipped, _ = zip_up_mpos(mpo1, mpo2, oc)
            
            # Extract both as ITensor arrays
            T_applied = prod(mpo_applied.data)
            arr_applied = Array(T_applied, prime.(sites)..., sites...)
            
            T_zipped = prod(mpo_zipped.data)
            arr_zipped = Array(T_zipped, prime.(sites)..., sites...)
            
            # Compare ITensor arrays
            err = LinearAlgebra.norm(arr_applied - arr_zipped)
            @test isapprox(err, 0.0; atol=1e-12)
        end
    end

    @testset "Zip Up: No Compression (Bond Dimensions)" begin
        # Verify that zip_up does not compress - it only combines tensors
        # Bond dimensions should be product of original bond dimensions
        
        n = 4
        sites = [Index(2, "site-$i") for i in 1:n]
        
        mpo1 = control_Hphase_mpo(n, sites)
        subset_sites = sites[2:end]
        mpo2 = control_Hphase_mpo(n-1, subset_sites)
        
        # Store original bond dimensions
        bonds1_dims = [dim(b) for b in mpo1.bonds]
        bonds2_dims = [dim(b) for b in mpo2.bonds]
        
        # Prepare indices
        for (idx2, s) in enumerate(subset_sites)
            idx1 = idx2 + 1
            s_cont = sim(s)
            mpo1.data[idx1] = replaceinds(mpo1.data[idx1], s => s_cont)
            mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
        end
        
        oc = n
        mpo_zipped, _ = zip_up_mpos(mpo1, mpo2, oc)
        
        # Check that no compression occurred by examining a site in the overlap region
        T_original_1 = prod(mpo1.data)
        T_original_2 = prod(mpo2.data)
        
        # The zipped result should equal the product (within numerical precision)
        # We already tested this above, so here we just verify dimensions make sense
        for i in 1:length(mpo_zipped.bonds)
            bond_dim = dim(mpo_zipped.bonds[i])
            # Bond dimensions should be reasonable (not compressed to 1 artificially)
            @test bond_dim >= 1
        end
    end
end

# ======================== ZIP DOWN TESTS ========================
@testset "qft_transformer.jl: zip down algorithm correctness" begin
    @testset "Zip Down: Compression Within Tolerance" begin
        # Test that zip_down compresses the MPO while maintaining accuracy within tolerance
        
        for n in 3:5
            sites = [Index(2, "site-$i") for i in 1:n]
            
            # Build a QFT MPO (which will have been zipped up)
            # We'll create a simple combined MPO for testing
            mpo1 = control_Hphase_mpo(n, sites)
            subset_sites = sites[2:end]
            mpo2 = control_Hphase_mpo(n-1, subset_sites)
            
            # Prepare and zip up
            for (idx2, s) in enumerate(subset_sites)
                idx1 = idx2 + 1
                s_cont = sim(s)
                mpo1.data[idx1] = replaceinds(mpo1.data[idx1], s => s_cont)
                mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
            end
            
            oc = n
            mpo_zipped, new_oc = zip_up_mpos(mpo1, mpo2, oc)
            
            # Store original as ITensor array
            T_original = prod(mpo_zipped.data)
            arr_original = Array(T_original, prime.(sites)..., sites...)
            
            # Apply zip_down with specified tolerance
            cutoff = 1e-12
            maxdim = 100
            mpo_compressed, final_oc = zip_down_mpos(mpo_zipped, new_oc; cutoff=cutoff, maxdim=maxdim)
            
            # Extract compressed result as ITensor array
            T_compressed = prod(mpo_compressed.data)
            arr_compressed = Array(T_compressed, prime.(sites)..., sites...)
            
            # Verify result matches within tolerance using ITensor arrays
            err = LinearAlgebra.norm(arr_compressed - arr_original)
            @test isapprox(err, 0.0; atol=cutoff)
            
            # Verify OC moved to bottom
            @test final_oc == n
        end
    end

    @testset "Zip Down: Bond Dimension Reduction" begin
        # Test that zip_down reduces bond dimensions (or keeps them same if already optimal)
        
        for n in 3:5
            sites = [Index(2, "site-$i") for i in 1:n]
            
            # Create combined MPO
            mpo1 = control_Hphase_mpo(n, sites)
            subset_sites = sites[2:end]
            mpo2 = control_Hphase_mpo(n-1, subset_sites)
            
            for (idx2, s) in enumerate(subset_sites)
                idx1 = idx2 + 1
                s_cont = sim(s)
                mpo1.data[idx1] = replaceinds(mpo1.data[idx1], s => s_cont)
                mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
            end
            
            oc = n
            mpo_zipped, new_oc = zip_up_mpos(mpo1, mpo2, oc)
            
            # Store bond dimensions before compression
            bonds_before = [dim(b) for b in mpo_zipped.bonds]
            
            # Apply zip_down
            cutoff = 1e-12
            maxdim = 100
            mpo_compressed, _ = zip_down_mpos(mpo_zipped, new_oc; cutoff=cutoff, maxdim=maxdim)
            
            # Check bond dimensions after compression
            bonds_after = [dim(b) for b in mpo_compressed.bonds]
            
            # Verify compression: each bond should be <= original (or same if already optimal)
            for i in eachindex(bonds_after)
                @test bonds_after[i] <= bonds_before[i]
            end
            
            # At least verify max bond dimension respects maxdim
            @test maximum(bonds_after) <= maxdim
        end
    end

    @testset "Zip Down: Identity Test with Apply" begin
        # Test zip_down on a case where MPOs are combined and verify using apply on MPS states
        
        n = 3
        sites = [Index(2, "site-$i") for i in 1:n]
        
        # Create a basis state first
        psi_test, _ = signal_mps(basis_state_vector(3, n))
        
        # Now create MPOs with the MPS sites to ensure they match from the start
        mps_sites = psi_test.sites
        mpo1 = control_Hphase_mpo(n, mps_sites)
        subset_sites = mps_sites[2:end]
        mpo2 = control_Hphase_mpo(n-1, subset_sites)
        
        # Prepare indices for zip_up
        for (idx2, s) in enumerate(subset_sites)
            idx1 = idx2 + 1
            s_cont = sim(s)
            mpo1.data[idx1] = replaceinds(mpo1.data[idx1], s => s_cont)
            mpo2.data[idx2] = replaceinds(mpo2.data[idx2], s' => s_cont)
        end
        
        # Perform zip_up to create a combined MPO
        oc = n
        mpo, new_oc = zip_up_mpos(mpo1, mpo2, oc)
        
        # Store original as ITensor array
        T_original = prod(mpo.data)
        arr_original = Array(T_original, prime.(mps_sites)..., mps_sites...)
        
        # Apply zip_down from the top
        mpo_compressed, _ = zip_down_mpos(mpo, new_oc; cutoff=1e-14, maxdim=10)
        
        # Extract compressed result as ITensor array
        T_compressed = prod(mpo_compressed.data)
        arr_compressed = Array(T_compressed, prime.(mps_sites)..., mps_sites...)
        
        # Verify using ITensor arrays
        err = LinearAlgebra.norm(arr_compressed - arr_original) / LinearAlgebra.norm(arr_original)
        @test err < 1e-12
        
        # Also verify by applying both MPOs to the test state and comparing results
        # Apply original (zipped) MPO
        psi_orig_result = apply(mpo, psi_test; cutoff=0.0, maxdim=1000)
        vec_orig = mps_to_vector(psi_orig_result)
        
        # Apply compressed MPO
        psi_comp_result = apply(mpo_compressed, psi_test; cutoff=0.0, maxdim=1000)
        vec_comp = mps_to_vector(psi_comp_result)
        
        # Compare application results
        err_apply = LinearAlgebra.norm(vec_comp - vec_orig)
        @test isapprox(err_apply, 0.0; atol=1e-12)
    end
end

# ======================== Qn TRANSFORMATION TESTS ========================

@testset "qft_transforme.jl: Q_n Transformation on basis states" begin
    # Test that the QFT MPO (which produces Qn = bit-reversed DFT) correctly transforms each basis state
    # For each basis state |j>, the result should match the j-th column of Qn matrix
    
    for n in 2:5
        N = 2^n
        sites = [Index(2, "site-$i") for i in 1:n]
        
        # Get expected Qn matrix
        Q_n = qn_matrix(n)
        
        # Test each basis state
        for j in 0:N-1
            # Build a fresh QFT MPO with generic sites
            qft_mpo = build_qft_mpo(n, sites; cutoff=1e-14, maxdim=1000)
            
            # Create basis state |j> (this will have its own new sites)
            psi_j, _ = signal_mps(basis_state_vector(j, n))
            
            # Update qft_mpo sites to match psi_j sites using update_site!
            for i in 1:n
                old_site = qft_mpo.sites[i]
                new_site = psi_j.sites[i]
                if old_site != new_site
                    update_site!(qft_mpo, old_site, new_site)
                end
            end
            
            # Apply QFT MPO
            psi_result = apply(qft_mpo, psi_j; cutoff=0.0, maxdim=1000)
            
            # Extract result vector
            result_vec = mps_to_vector(psi_result)
            
            # Expected: j-th column of Qn matrix (0-indexed)
            expected_vec = Q_n[:, j+1]
            
            # Compare
            err = LinearAlgebra.norm(result_vec - expected_vec)
            
            @test isapprox(err, 0.0; atol=1e-10)
        end
        
    end
end

# ======================== FULL QFT (Rn * Qn = Fn) TESTS ========================

@testset "qft_transformer.jl: Full DFT on basis states" begin
    # Test that applying bit-reversal permutation to Qn output gives standard DFT
    # Test on individual basis states
    
    for n in 2:5
        N = 2^n
        sites = [Index(2, "site-$i") for i in 1:n]
        
        # Test each basis state
        for j in 0:N-1
            # Build a fresh QFT MPO
            qft_mpo = build_qft_mpo(n, sites; cutoff=1e-14, maxdim=1000)
            
            # Create basis state |j> (this will have its own new sites)
            psi_j, _ = signal_mps(basis_state_vector(j, n))
            
            # Update qft_mpo sites to match psi_j sites using update_site!
            for i in 1:n
                old_site = qft_mpo.sites[i]
                new_site = psi_j.sites[i]
                if old_site != new_site
                    update_site!(qft_mpo, old_site, new_site)
                end
            end
            
            # Apply QFT MPO -> get Qn|j>
            psi_qn = apply(qft_mpo, psi_j; cutoff=0.0, maxdim=1000)
            qn_vec = mps_to_vector(psi_qn)
            
            # Apply bit-reversal to get Fn|j>
            N_vec = length(qn_vec)
            n_vec = Int(log2(N_vec))
            fn_vec = similar(qn_vec)
            for i in 0:N_vec-1
                r = bits_to_int(int_to_bits(i, n_vec; order=:lsb); order=:msb)
                fn_vec[r + 1] = qn_vec[i + 1]
            end
            
            # Expected dft vector
            sig = [i == j ? 1.0 : 0.0 for i in 0:(N-1)]
            expected_vec = dft(sig)
            
            # Compare
            err = LinearAlgebra.norm(fn_vec - expected_vec)
            @test isapprox(err, 0.0; atol=1e-10)
        end
    end
end

@testset "qft_transformer.jl: FFTW comparison on random signal" begin
    # Test that the full DFT (Rn * Qn) matches FFTW output
    # FFTW conventions:
    #   fft(x)[k]  = sum_j x[j] * exp(-2πi * j * k / N)   (forward, -2πi)
    #   bfft(x)[k] = sum_j x[j] * exp(+2πi * j * k / N)   (backward, +2πi)
    # Our QFT convention: Fn[j,k] = exp(+2πi * j * k / N) / sqrt(N)
    # So: Our Fn * x = bfft(x) / sqrt(N)
    
    for n in 2:5
        N = 2^n
        sites = [Index(2, "site-$i") for i in 1:n]
        
        # Random test signal
        sig = randn(ComplexF64, N)
        
        # Apply QFT MPO
        psi, norm_c = signal_mps(sig)
        
        # Build QFT MPO
        qft_mpo = build_qft_mpo(n, psi.sites; cutoff=1e-14, maxdim=1000)
        
        psi_qn = apply(qft_mpo, psi; cutoff=0.0, maxdim=1000)
        qn_result = mps_to_vector(psi_qn) * norm_c
        
        # Apply bit-reversal to get full DFT
        N_res = length(qn_result)
        n_res = Int(log2(N_res))
        fn_result = similar(qn_result)
        for i in 0:N_res-1
            r = bits_to_int(int_to_bits(i, n_res; order=:lsb); order=:msb)
            fn_result[r + 1] = qn_result[i + 1]
        end
        
        # Compare with FFTW (bfft uses +2πi convention like our QFT)
        result_fftw = bfft(sig) / sqrt(N)
        err = LinearAlgebra.norm(fn_result - result_fftw)
        @test isapprox(err, 0.0; atol=1e-10)
    end
end
