import QILaplace.QFTGates: I, H, P, Π, control_Hphase_mpo

# Elementary circuit gates testing
@testset "qft_gates.jl: Elementary Quantum circuit gates correctness" begin
    # Test element-by-element correctness of I, H, P, Π gates
    site = Index(2, "site")

    I_gate = I(site)
    for (i, j) in Iterators.product(1:2, 1:2)
        expected = i == j ? 1.0 : 0.0
        @test isapprox(I_gate[i, j], expected; atol=1e-12, rtol=1e-12)
    end
    
    H_gate = H(site)
    H_expected = (1/sqrt(2)) * [1 1; 1 -1]
    for (i, j) in Iterators.product(1:2, 1:2)
        @test isapprox(H_gate[i, j], H_expected[i, j]; atol=1e-12, rtol=1e-12)
    end

    θs = [0.0, π/2, π, 3π/2]
    for θ in θs
        P_gate = P(θ, site)
        P_expected = [1 0; 0 exp(im * θ)]
        for (i, j) in Iterators.product(1:2, 1:2)
            @test isapprox(P_gate[i, j], P_expected[i, j]; atol=1e-12, rtol=1e-12)
        end
    end

    for i in 0:1
        Π_gate = Π(i, site)
        for (m, n) in Iterators.product(1:2, 1:2)
            expected = (m == n && m == (i + 1)) ? 1.0 : 0.0
            @test isapprox(Π_gate[m, n], expected; atol=1e-12, rtol=1e-12)
        end
    end
end

# test the controlled phase gate MPO by comparing the matrix multiplication of the control phase gates with the dense representation of the MPO
@testset "qft_gates.jl: control_Hphase_mpo matches analytical basis action" begin
    for n in 2:4
        sites = [Index(2, "site-$i") for i in 1:n]
        W = control_Hphase_mpo(n, sites)

        for b in 0:((1 << n) - 1)
            # Create input basis state using signal_mps
            x = [i == (b+1) ? 1.0 : 0.0 for i in 1:(1 << n)]
            ψ_in, _ = signal_mps(x)
            
            # Replace signal_mps sites with MPO sites
            for i in 1:n
                update_site!(ψ_in, ψ_in.sites[i], sites[i])
            end

            # Apply MPO and extract dense output vector
            ψ_out = apply(W, ψ_in)
            T_out = to_dense_mps(ψ_out)
            A_out = Array(T_out, sites...)
            v_out = ComplexF64.(vec(A_out))

            # Analytical expected output
            site_bits_input = _int_to_bit(b, n)
            
            # Compute total phase for the |site[1]=1⟩ component
            phase_for_1 = 1.0 + 0.0im
            for l in 2:n
                if site_bits_input[l] == 1  # apply phase when site[l] is |1⟩
                    θ = 2π / 2.0^l
                    phase_for_1 *= exp(im * θ)
                end
            end

            # The expected output vector
            # Amplitude on |site[1]=0, site[2..n] unchanged⟩
            a0 = (site_bits_input[1] == 0 ? 1 : 1) / sqrt(2)
            # Amplitude on |site[1]=1, site[2..n] unchanged⟩ with accumulated phase
            a1 = (site_bits_input[1] == 0 ? 1 : -1) / sqrt(2) * phase_for_1

            # Output states: |site[1]=0, others unchanged> and |site[1]=1, others unchanged>
            idx_rest = 0
            for i in 2:n
                idx_rest += site_bits_input[i] << (i - 1)
            end
            idx0 = idx_rest      # site[1] = 0, others from input
            idx1 = idx_rest + 1  # site[1] = 1, others from input

            v_exp = zeros(ComplexF64, 1 << n)
            v_exp[idx0 + 1] = a0
            v_exp[idx1 + 1] = a1

            @test isapprox(v_out, v_exp; atol=1e-10)
        end
    end
end
