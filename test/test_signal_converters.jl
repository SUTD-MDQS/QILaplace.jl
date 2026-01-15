using QILaplace.SignalConverters: _array_to_tensor, _tensor_to_mps_svd, _tensor_to_mps_rsvd, _tensor_to_mps

# ==================== Helper: Test _array_to_tensor ====================
@testset "SignalConverters.jl: _array_to_tensor correctness" begin
    # Test with power-of-2 length signal
    x = [1.0, 2.0, 3.0, 4.0]
    ψ, normalisation = _array_to_tensor(x)
    @test length(inds(ψ)) == 2
    @test all(dim.(inds(ψ)) .== 2)
    @test isapprox(normalisation, sqrt(1.0 + 4.0 + 9.0 + 16.0); rtol=1e-12)
    
    # Verify normalized tensor
    arr = Array(ψ, inds(ψ)...)
    @test isapprox(sum(abs2, arr), 1.0; rtol=1e-12)
    
    # Test with non-power-of-2 (should warn and pad)
    x_non_pow2 = [1.0, 2.0, 3.0]
    ψ_padded, norm_padded = _array_to_tensor(x_non_pow2)
    @test length(inds(ψ_padded)) == 2  # rounds up to 4 = 2^2
    @test isapprox(norm_padded, sqrt(1.0 + 4.0 + 9.0); rtol=1e-12)
end

# ==================== Test _tensor_to_mps_svd ====================
@testset "SignalConverters.jl: _tensor_to_mps_svd reconstruction accuracy" begin
    n = 5
    x = randn(2^n)
    ψ, _ = _array_to_tensor(x)
    
    # Convert to MPS using SVD
    mps_svd = _tensor_to_mps_svd(ψ; cutoff=1e-14, maxdim=typemax(Int))
    
    # Contract MPS back to tensor
    T_reconstructed = contract_chain(mps_svd.data)
    arr_orig = Array(ψ, inds(ψ)...)
    arr_recon = Array(T_reconstructed, inds(ψ)...)
    
    # Check reconstruction accuracy
    @test isapprox(arr_orig, arr_recon; atol=1e-12, rtol=1e-12)
    
    # Verify MPS structure
    @test length(mps_svd.data) == n
    @test length(mps_svd.sites) == n
    @test length(mps_svd.bonds) == n - 1
    
    # Test single-site edge case
    x_single = [0.6, 0.8]
    ψ_single, _ = _array_to_tensor(x_single)
    mps_single = _tensor_to_mps_svd(ψ_single)
    @test length(mps_single.data) == 1
    @test length(mps_single.bonds) == 0
end

# ==================== Test _tensor_to_mps_rsvd ====================
@testset "SignalConverters.jl: _tensor_to_mps_rsvd reconstruction accuracy" begin
    n = 5
    x = randn(2^n)
    ψ, _ = _array_to_tensor(x)
    
    # Convert to MPS using RSVD
    mps_rsvd = _tensor_to_mps_rsvd(ψ; cutoff=1e-12, maxdim=128, k=32, p=5, q=2)
    
    # Contract MPS back to tensor
    T_reconstructed = contract_chain(mps_rsvd.data)
    arr_orig = Array(ψ, inds(ψ)...)
    arr_recon = Array(T_reconstructed, inds(ψ)...)
    
    # Check reconstruction accuracy (RSVD is approximate, so slightly looser tolerance)
    @test isapprox(arr_orig, arr_recon; atol=1e-10, rtol=1e-10)
    
    # Verify MPS structure
    @test length(mps_rsvd.data) == n
    @test length(mps_rsvd.sites) == n
    @test length(mps_rsvd.bonds) == n - 1
    
    # Test single-site edge case
    x_single = [0.6, 0.8]
    ψ_single, _ = _array_to_tensor(x_single)
    mps_single = _tensor_to_mps_rsvd(ψ_single)
    @test length(mps_single.data) == 1
    @test length(mps_single.bonds) == 0
end

# ==================== Compare SVD vs RSVD ====================
@testset "SignalConverters.jl: SVD vs RSVD comparison" begin
    n = 6
    x = randn(2^n)
    ψ, _ = _array_to_tensor(x)

    mps_svd = _tensor_to_mps_svd(ψ; cutoff=1e-12, maxdim=128)
    mps_rsvd = _tensor_to_mps_rsvd(ψ; cutoff=1e-12, maxdim=128, k=64, p=5, q=2)

    # Reconstruct full tensors
    T_svd = contract_chain(mps_svd.data)
    T_rsvd = contract_chain(mps_rsvd.data)

    arr_orig = Array(ψ, inds(ψ)...)
    arr_svd = Array(T_svd, inds(ψ)...)
    arr_rsvd = Array(T_rsvd, inds(ψ)...)

    # Check approximate equality to original
    @test isapprox(arr_orig, arr_svd; atol=1e-10, rtol=1e-10)
    @test isapprox(arr_orig, arr_rsvd; atol=1e-9, rtol=1e-9)
    
    # Check SVD and RSVD are close to each other
    @test isapprox(arr_svd, arr_rsvd; atol=1e-8, rtol=1e-8)

    # Bond counts
    @test length(mps_svd.bonds) == n - 1
    @test length(mps_rsvd.bonds) == n - 1
end

# ==================== Test signal_mps (external function) ====================
@testset "SignalConverters.jl: signal_mps with various methods" begin
    x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]  # n=3
    
    # Test SVD method
    ψ_svd, norm_svd = signal_mps(x; method=:svd, cutoff=1e-14)
    @test isapprox(norm_svd, sqrt(sum(abs2, x)); rtol=1e-12)
    
    # Reconstruct and verify coefficient-by-coefficient
    for i in 0:7
        expected = x[i+1] / norm_svd
        obtained = coefficient(ψ_svd, i)
        @test isapprox(obtained, expected; atol=1e-12)
    end
    
    # Test RSVD method
    ψ_rsvd, norm_rsvd = signal_mps(x; method=:rsvd, cutoff=1e-12, maxdim=128)
    @test isapprox(norm_rsvd, norm_svd; rtol=1e-12)
    
    # Compare RSVD to SVD via coefficients
    for i in 0:7
        svd_coeff = coefficient(ψ_svd, i)
        rsvd_coeff = coefficient(ψ_rsvd, i)
        @test isapprox(svd_coeff, rsvd_coeff; atol=1e-10, rtol=1e-10)
    end
    
    # Test error handling
    @test_throws ArgumentError signal_mps(x; method=:invalid_method)
end

# ==================== Coefficient extraction: various input types ====================
@testset "SignalConverters.jl: coefficient extraction with all input types" begin
    # Create a simple signal with known values
    x = Float64[1, 2, 3, 4, 5, 6, 7, 8]  # n=3, indices 000 to 111
    ψ, normalisation = signal_mps(x; method=:svd, cutoff=0.0, maxdim=typemax(Int))
    
    # Test 1: Vector of integers (zero-based)
    for i in 0:7
        bits = [(i >> (2-j)) & 1 for j in 0:2]  # Extract 3 bits
        expected = x[i+1] / normalisation
        obtained = coefficient(ψ, bits)
        @test isapprox(obtained, expected; atol=1e-12)
    end
    
    # Test 2: Direct integer input (big-endian)
    @test isapprox(coefficient(ψ, 0), x[1] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, 1), x[2] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, 2), x[3] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, 5), x[6] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, 7), x[8] / normalisation; atol=1e-12)
    
    # Test 3: Bitstring (contiguous)
    @test isapprox(coefficient(ψ, "000"), x[1] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, "001"), x[2] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, "010"), x[3] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, "101"), x[6] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, "111"), x[8] / normalisation; atol=1e-12)
    
    # Test 4: Bitstring with separators
    @test isapprox(coefficient(ψ, "[0,0,0]"), x[1] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, "1, 0, 1"), x[6] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, "[1,1,1]"), x[8] / normalisation; atol=1e-12)
    
    # Test 5: getindex syntax (ψ[...])
    @test isapprox(ψ[0, 0, 0], x[1] / normalisation; atol=1e-12)
    @test isapprox(ψ[0, 0, 1], x[2] / normalisation; atol=1e-12)
    @test isapprox(ψ[1, 0, 1], x[6] / normalisation; atol=1e-12)
    @test isapprox(ψ[1, 1, 1], x[8] / normalisation; atol=1e-12)
    
    # Test 6: Tuple input
    @test isapprox(coefficient(ψ, (0, 0, 0)), x[1] / normalisation; atol=1e-12)
    @test isapprox(coefficient(ψ, (1, 0, 1)), x[6] / normalisation; atol=1e-12)
    
    # Test error cases
    @test_throws ArgumentError coefficient(ψ, [0, 0])  # wrong length
    @test_throws ArgumentError coefficient(ψ, [2, 0, 0])  # out of range
    @test_throws ArgumentError coefficient(ψ, 0b1000)  # integer too large (needs 4 bits)
end

# ==================== RSVD coefficient validation ====================
@testset "SignalConverters.jl: RSVD coefficient accuracy" begin
    x = Float64[1, 2, 3, 4, 5, 6, 7, 8]
    ψ_rsvd, normalisation = signal_mps(x; method=:rsvd, cutoff=1e-12, maxdim=128)
    
    # Test with various input types
    @test isapprox(ψ_rsvd[0, 0, 0], x[1] / normalisation; atol=1e-10)
    @test isapprox(coefficient(ψ_rsvd, 0b010), x[3] / normalisation; atol=1e-10)
    @test isapprox(coefficient(ψ_rsvd, "101"), x[6] / normalisation; atol=1e-10)
    @test isapprox(coefficient(ψ_rsvd, [1,1,1]), x[8] / normalisation; atol=1e-10)
end
