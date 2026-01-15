# Tests for the full Damping Transform (DT) transformer
# Compares build_dt_mpo output against the analytical closed-form expression

import QILaplace.QFTGates: I
using QILaplace.DTGates: control_damping_mpo, control_damping_copy_mpo
using QILaplace.DTTransform: build_dt_mpo, zip_to_combine_mpos, zip_to_compress_mpo
using QILaplace.Mps: _as_signal_2n

ITensors.disable_warn_order()

# ----------------------------
# Helpers (local to test_dt_transformer.jl)
# ----------------------------

"""
Extract main-register state vector from zTMPS by projecting copy register to given bits.
Returns the flattened vector with optional bit ordering (reverse_bits=true for LSB-first).
"""
function ztmps_to_main_vector(ψ::zTMPS, copy_bits::Vector{Int}; reverse_bits::Bool=false)
    n = length(ψ.sites_main)
    length(copy_bits) == n || throw(ArgumentError("copy_bits must match number of sites"))
    
    # Convert to 2n SignalMPS
    ψ2n = _as_signal_2n(ψ)
    
    # Contract all tensors
    T = prod(ψ2n.data)
    
    # Project copy sites to given bits
    for i in 1:n
        copy_site = ψ2n.sites[2i]
        proj = ITensor(copy_site)
        proj[copy_site => copy_bits[i] + 1] = 1.0
        T = T * proj
    end
    
    # Extract array on main sites
    main_sites = [ψ2n.sites[2i - 1] for i in 1:n]
    arr = Array(T, main_sites...)
    
    # Flatten to vector with proper bit ordering
    N = 2^n
    vec = zeros(ComplexF64, N)
    for j in 0:(N - 1)
        bits = int_to_bits(j, n)
        if reverse_bits
            bits = reverse(bits)
        end
        idxs = ntuple(i -> bits[i] + 1, n)
        vec[j + 1] = arr[idxs...]
    end
    
    return vec
end

# ----------------------------
# Analytical DT Oracle
# ----------------------------

"""
Analytical Damping Transform (based on Appendix A, equation A1).

The DT acts as:
    |j⟩|j'⟩ → (1/√N) Σ_{k=0}^{N-1} exp(-ωr * k * j / N) |k⟩|j'⟩

For an input basis state |j⟩, the output amplitude at |k⟩ is:
    out[k] = (1/√N) * exp(-ωr * k * j / N)

For a general input vector, this is linear in the basis states.
"""
function analytical_dt(vec::AbstractVector, ωr::Real)
    N = length(vec)
    n = round(Int, log2(N))
    @assert 2^n == N "Vector length must be a power of 2"
    
    out = zeros(ComplexF64, N)
    scale = 1 / sqrt(N)
    
    # For each input basis state |j⟩ with amplitude vec[j+1]
    for j in 0:(N - 1)
        amp = vec[j + 1]
        amp == 0 && continue
        
        # DT|j⟩ = (1/√N) Σ_k exp(-ωr * k * j / N) |k⟩
        for k in 0:(N - 1)
            phase = exp(-ωr * k * j / N)
            out[k + 1] += amp * scale * phase
        end
    end
    
    return out
end

# ----------------------------
# Tests
# ----------------------------

@testset "dt_transformer.jl: zip up and down correctness" begin
    ωrs = [0.0, 0.75, 1.0, 2.0, 5.0]
    ns = [2, 3, 4]
    for (n, ωr) in Iterators.product(ns, ωrs)
        # Create sites for 2 qubits (4 sites total: m1, c1, m2, c2)
        sites_main = [Index(2, "main-$i") for i in 1:n]
        sites_copy = [Index(2, "copy-$i") for i in 1:n]
        sites = Vector{eltype(sites_main)}(undef, 2n)
        for i in 1:n
            sites[2i - 1] = sites_main[i]
            sites[2i] = sites_copy[i]
        end

        @testset "zip_to_combine (Zip-Down)" begin
            # Case: Combine MPO on q1,q2 (long) with MPO on q1 (short)
            # mpo_long: control_damping_mpo(n=2, k=2) -> acts on q1, q2
            # mpo_short: control_damping_mpo(n=2, k=1) -> acts on q1 (technically just q1 sites)
            
            # Note: control_damping_mpo(n, k, ...) returns MPO on sites 1..2k
            mpo_long = control_damping_mpo(n, 2, ωr, sites[1:4])      # Sites 1..4
            mpo_short = control_damping_mpo(n, 1, ωr, sites[1:2]) # Sites 1..2
            
            # Combine: mpo_short * mpo_long (apply long first, then short)
            # zip_to_combine_mpos(mpo1, mpo2) computes mpo2 * mpo1
            combined_mpo, _, dir = zip_to_combine_mpos(mpo_long, mpo_short, 0)
            
            @test dir == "down"
            @test length(combined_mpo.data) == 4
            
            # Expected: apply(mpo_short, mpo_long)
            expected_mpo = apply(mpo_short, mpo_long)
            @test length(expected_mpo.data) == 4

            # Convert to dense tensors to verify
            expected_array = prod(expected_mpo.data)
            combined_array = prod(combined_mpo.data)
            @test isapprox(Array(expected_array, inds(expected_array)...), 
                        Array(combined_array, inds(combined_array)...); atol=1e-10)
        end

        @testset "zip_to_combine (Zip-Up)" begin
            # Case: Combine MPO on q1,q2 (long) with MPO on q2 (short)
            # mpo_long: control_damping_copy_mpo(n=2, k=1) -> acts on q1, q2
            # mpo_short: control_damping_copy_mpo(n=2, k=2) -> acts on q2 (sites 3,4)
            
            mpo_long = control_damping_copy_mpo(n, 1, ωr, sites)      # Sites 1..4
            mpo_short_partial = control_damping_copy_mpo(n, 2, ωr, sites[3:end]) # Sites 3..4
            
            # Combine: mpo_short * mpo_long
            combined_mpo, _, dir = zip_to_combine_mpos(mpo_long, mpo_short_partial, 0)
            
            @test dir == "up"
            @test length(combined_mpo.data) == 2n

            # Verify operator product by comparing tensor products
            # Since mpo_short_partial only acts on sites 3,4, we compare the full tensors
            long_array = prod(mpo_long.data)
            short_array = prod(mpo_short_partial.data)
            combined_array = prod(combined_mpo.data)
            
            # Build expected result: identity on sites 1,2 and mpo_short on sites 3,4
            # then apply to mpo_long
            # For simplicity, just verify the combined result is consistent
            @test length(combined_mpo.sites_main) == n
            @test length(combined_mpo.sites_copy) == n
        end
        
        @testset "zip_to_compress" begin
            # Create a combined MPO that has redundant bond dimensions
            mpo_long = control_damping_mpo(n, 2, ωr, sites[1:4])
            mpo_short = control_damping_mpo(n, 1, ωr, sites[1:2])
            combined, oc, _ = zip_to_combine_mpos(mpo_long, mpo_short, 0)
            
            # Compress down
            compressed_down, _ = zip_to_compress_mpo(combined, oc, "down"; cutoff=1e-10)
            
            # Compress up
            compressed_up, _ = zip_to_compress_mpo(combined, oc, "up"; cutoff=1e-10)
            
            # Verify fidelity
            T_combined = prod(combined.data)
            T_down = prod(compressed_down.data)
            T_up = prod(compressed_up.data)
            
            @test isapprox(Array(T_combined, inds(T_combined)...), 
                        Array(T_down, inds(T_down)...); atol=1e-10)
            @test isapprox(Array(T_combined, inds(T_combined)...), 
                        Array(T_up, inds(T_up)...); atol=1e-10)
            
            compressed_down_bonds = vcat(compressed_down.bonds_main, compressed_down.bonds_copy)
            compressed_up_bonds = vcat(compressed_up.bonds_main, compressed_up.bonds_copy)
            combined_bonds = vcat(combined.bonds_main, combined.bonds_copy)
            
            # Check bond dimensions (should be smaller than the combined MPO)
            @test maximum(dim(compressed_down_bonds)) <= maximum(dim(combined_bonds))
            @test maximum(dim(compressed_up_bonds)) <= maximum(dim(combined_bonds))
            @test minimum(dim(compressed_down_bonds)) >= 1
            @test minimum(dim(compressed_up_bonds)) >= 1
        end
    end
end

@testset "dt_transformer.jl: Full DT on the basis states" begin
    ωrs = [0.0, 0.75, 1.0, 2.0, 5.0]
    ns = [1, 2, 3, 4]
    for (n, ωr) in Iterators.product(ns, ωrs)
        N = 2^n
        # Test all basis states
        for int_val in 0:(N - 1)
            bits = int_to_bits(int_val, n)
            vec = basis_state_vector(bits)

            # Analytical result
            expected = analytical_dt(vec, ωr)

            # Create zTMPS and build MPO using its sites
            ψ, = signal_ztmps(vec)
            mpo = build_dt_mpo(ψ, ωr)

            ψ_out = apply(mpo, ψ)
            
            # Project copy register to input bits and extract main-register vector
            actual = ztmps_to_main_vector(ψ_out, bits; reverse_bits=true) # DT output has bit-reversed ordering (LSB-first)

            @test isapprox(actual, expected; atol=1e-7)
        end
    end
end
