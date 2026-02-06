import QILaplace.DTGates: I, dampedH, R, Π, control_damping_mpo, control_damping_copy_mpo
import QILaplace.Mps: _as_signal_2n

# Elementary circuit gates testing
@testset "dt_gates.jl: Elementary non-unitary circuit gates correctness" begin
    # Test element-by-element correctness of I, dampedH, R, Π gates
    site = Index(2, "site")

    ωrs = [0.0, 0.25, 0.5, 1.0, 1.1]
    for ωr in ωrs
        dampedH_gate = dampedH(ωr, site)
        dampedH_expected = 1/sqrt(2) * [1.0 1.0; 1.0 exp(-ωr/2)]
        for (i, j) in Iterators.product(1:2, 1:2)
            @test isapprox(
                dampedH_gate[i, j], dampedH_expected[i, j]; atol=1e-12, rtol=1e-12
            )
        end

        R_gate = R(ωr, site)
        R_expected = [1.0 0.0; 0.0 exp(-ωr)]
        for (i, j) in Iterators.product(1:2, 1:2)
            @test isapprox(R_gate[i, j], R_expected[i, j]; atol=1e-12, rtol=1e-12)
        end
    end
end

################################# HELPER FUNCTIONS #####################################

"""Extract main and copy bits from interleaved bit vector."""
function split_main_copy_bits(bits_2n::Vector{Int})
    n = length(bits_2n) ÷ 2
    main_bits = [bits_2n[2i - 1] for i in 1:n]
    copy_bits = [bits_2n[2i] for i in 1:n]
    return main_bits, copy_bits
end

################################# CONTROL_DAMPING_MPO TESTS #####################################

# NOTE ON BIT ORDERING CONVENTION:
# signal_ztmps(x) creates a zTMPS from flat vector x where x[b+1] corresponds to basis state |b⟩
# The zTMPS interleaves main and copy sites: [main[1], copy[1], main[2], copy[2], ...]
# When extracted with Array(T, sites...), site[1] is fastest-varying (LSB in flat index).
#
# For the 2n-site interleaved representation with sites[2i-1]=main[i], sites[2i]=copy[i]:
# flat index k corresponds to: k = Σ site[i] * 2^(i-1)
# where site[2i-1] is main[i] and site[2i] is copy[i]

import QILaplace.ApplyMPO: apply, _as_single_site_mpo

@testset "dt_gates.jl: control_damping_mpo matches analytical basis action" begin
    ωrs = [0.0, 0.5, 1.0, 2.0, 5.0] # damping rates
    ks = [1, 2, 3]  # number of qubits in the controlled window

    for (ωr, k) in Iterators.product(ωrs, ks)
        # Create sites for k paired qubits
        sites_main = [Index(2, "main-$i") for i in 1:k]
        sites_copy = [Index(2, "copy-$i") for i in 1:k]
        sites = vcat([[sites_main[i], sites_copy[i]] for i in 1:k]...)

        # Create MPO
        W = control_damping_mpo(k, k, ωr, sites)

        # Test all basis states
        for b in 0:((1 << k) - 1)
            # Create input basis state using signal_ztmps
            x = [i == (b+1) ? 1.0 : 0.0 for i in 1:(1 << k)]
            ψ_in, _ = signal_ztmps(x)

            # Replace signal_ztmps sites with MPO sites
            for i in 1:k
                update_site!(ψ_in, ψ_in.sites_main[i], sites_main[i])
                update_site!(ψ_in, ψ_in.sites_copy[i], sites_copy[i])
            end

            # Apply MPO and extract dense output vector
            ψ_out = apply(W, ψ_in)
            v_out = mps_to_vector(ψ_out)

            # Analytical expected output                
            bits = int_to_bits(b, k)  # bits[1] is MSB (qubit 1), bits[k] is LSB (qubit k)

            # The control is at qubit k (LSB)
            # The dampedH gate at site k creates a superposition:
            #   |0⟩ → (1/√2)|0⟩ + (1/√2)|1⟩  (when control=0, via Π0)
            #   |1⟩ → (1/√2)|0⟩ + (1/√2×exp(-ωr/2))|1⟩  (when control=1, via Π1)

            v_exp = zeros(ComplexF64, 1 << (2k))

            control_bit = bits[k]  # qubit k is the control (LSB)
            target_bits = bits[1:(k - 1)]  # qubits 1..k-1 are targets (MSBs)

            # The MPO only modifies MAIN qubits, COPY qubits stay unchanged from input
            copy_bits = bits  # Copy bits unchanged from input

            if control_bit == 0
                for hadamard_out in 0:1
                    # Π0 branch: Hadamard without damping on |0⟩ output
                    amp = 1.0 / sqrt(2.0)

                    # Targets: I gates, so main bits stay the same
                    main_out_bits = copy(target_bits)
                    push!(main_out_bits, hadamard_out)  # Append control qubit output

                    # Convert to zTMPS index (interleaved main/copy)
                    idx_out = 0
                    for l in 1:k
                        idx_out += main_out_bits[l] * (1 << (2*(l-1)))      # main[l] bit at position 2*(l-1)
                        idx_out += copy_bits[l] * (1 << (2*(l-1) + 1))      # copy[l] bit at position 2*(l-1)+1
                    end

                    v_exp[idx_out + 1] = amp
                end
            else
                # Π1 branch: Hadamard with damping on |1⟩ output
                for hadamard_out in 0:1
                    # Hadamard amplitude
                    if hadamard_out == 0
                        amp = 1.0 / sqrt(2.0)
                    else
                        amp = 1.0 / sqrt(2.0) * exp(-ωr / 2.0)
                    end

                    # Targets: R gates apply based on target bit values
                    main_out_bits = copy(target_bits)
                    for j in 1:(k - 1)
                        if main_out_bits[j] == 1
                            # R gate damping factor for qubit j (target position j, control at k)
                            θ_j = ωr * 2.0^(j - k - 1)
                            amp *= exp(-θ_j)
                        end
                    end
                    push!(main_out_bits, hadamard_out)  # Append control qubit output

                    # Convert to zTMPS index (interleaved main/copy)
                    idx_out = 0
                    for l in 1:k
                        idx_out += main_out_bits[l] * (1 << (2*(l-1)))      # main[l] bit at position 2*(l-1)
                        idx_out += copy_bits[l] * (1 << (2*(l-1) + 1))      # copy[l] bit at position 2*(l-1)+1
                    end

                    v_exp[idx_out + 1] = amp
                end
            end

            @test isapprox(v_out, v_exp; atol=1e-10)
        end
    end
end

################################# CONTROL_DAMPING_COPY_MPO TESTS #####################################

# control_damping_copy_mpo is a diagonal gate:
# - Control: copy[1] (copy of the k-th qubit in global terms)
# - Targets: main[2..L] where L = n - k + 1
# - Applies R(θ) gates diagonally: R(θ_j) on main[j] when copy[1]=|1⟩
# - No superposition is created (unlike control_damping_mpo with dampedH)

@testset "dt_gates.jl: control_damping_copy_mpo matches analytical basis action" begin
    ωrs = [0.0, 0.5, 1.0, 2.0, 5.0]

    for ωr in ωrs
        for n in 2:4
            for k in 1:(n - 1)  # k must be < n for non-trivial gate
                L = n - k + 1  # number of qubits in the window

                # Create sites for L paired qubits
                sites_main = [Index(2, "main-$i") for i in 1:L]
                sites_copy = [Index(2, "copy-$i") for i in 1:L]
                sites = vcat([[sites_main[i], sites_copy[i]] for i in 1:L]...)

                # Create MPO
                W = control_damping_copy_mpo(n, k, ωr, sites)

                # Iterate over L-qubit basis states
                for b in 0:((1 << L) - 1)
                    # Create input basis state using signal_ztmps
                    x = [i == (b+1) ? 1.0 : 0.0 for i in 1:(1 << L)]
                    ψ_in, _ = signal_ztmps(x)

                    # Replace signal_ztmps sites with MPO sites
                    for i in 1:L
                        update_site!(ψ_in, ψ_in.sites_main[i], sites_main[i])
                        update_site!(ψ_in, ψ_in.sites_copy[i], sites_copy[i])
                    end

                    # Apply MPO and extract dense output vector
                    ψ_out = apply(W, ψ_in)
                    v_out = mps_to_vector(ψ_out)

                    # Analytical expected output
                    # This is a DIAGONAL gate, so input basis state → scalar × same basis state
                    # Control is copy[1], targets are main[2..L]

                    bits = int_to_bits(b, L)  # bit[1] is MSB, bit[L] is LSB

                    # Input state is |bits⟩_main ⊗ |bits⟩_copy
                    # copy[1] corresponds to bit[1] (same as main[1] in input)

                    # Compute damping factor
                    damping = 1.0
                    if bits[1] == 1  # control is on
                        for j in 2:L
                            if bits[j] == 1  # target is |1⟩
                                θ_j = ωr * 2.0^(j - 2)
                                damping *= exp(-θ_j)
                            end
                        end
                    end

                    # Output is same basis state scaled by damping factor
                    idx_out = 0
                    for l in 1:L
                        idx_out += bits[l] * (1 << (2l - 2))  # main bit
                        idx_out += bits[l] * (1 << (2l - 1))  # copy bit
                    end

                    v_exp = zeros(ComplexF64, 1 << (2L))
                    v_exp[idx_out + 1] = damping

                    @test isapprox(v_out, v_exp; atol=1e-10)
                end
            end
        end
    end
end
