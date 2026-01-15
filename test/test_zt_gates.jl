import QILaplace.ZTGates: control_Hphase_ztmps_mpo
import QILaplace.QFTGates: H, P, I, Π
import QILaplace.Mps: _as_signal_2n

################################# CONTROL_HPHASE_ZTMPS_MPO TESTS #####################################

@testset "zt_gates.jl: control_Hphase_ztmps_mpo matches analytical basis action" begin
    for k in 1:3
        # Create sites for k paired qubits
        sites_main = [Index(2, "main-$i") for i in 1:k]
        sites_copy = [Index(2, "copy-$i") for i in 1:k]
        sites = vcat([[sites_main[i], sites_copy[i]] for i in 1:k]...)
        
        # Create MPO
        W = control_Hphase_ztmps_mpo(k, sites)
        
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
            
            # The control is at copy qubit k (LSB)
            #   control output |0⟩ → bond 1 → I gates on copy targets (no phase)
            #   control output |1⟩ → bond 2 → P gates on copy targets (with phase)
            
            v_exp = zeros(ComplexF64, 1 << (2k))
            
            control_bit = bits[k]  # qubit k is the control (LSB)
            target_bits = bits[1:k-1]  # qubits 1..k-1 are targets (MSBs)
            
            # The MPO only modifies COPY qubits, MAIN qubits stay unchanged from input
            main_bits = bits
            
            # Hadamard on control (copy qubit k)
            # For |0⟩: (1/√2)|0⟩ + (1/√2)|1⟩
            # For |1⟩: (1/√2)|0⟩ - (1/√2)|1⟩  (standard Hadamard)
            
            # Phase accumulation based on target INPUT bits
            phase_factor = ComplexF64(1.0)
            if control_bit == 1
                for j in 1:(k-1)
                    if target_bits[j] == 1
                        # Phase angle for copy site j
                        θ_j = -2π / 2.0^(k - j + 1)
                        phase_factor *= exp(im * θ_j)
                    end
                end
            end
            
            # Hadamard amplitudes with phase applied when control=1
            # Output |0⟩: 
            a0 = (control_bit == 0 ? 1.0 : 1.0) / sqrt(2.0) * phase_factor
            # Output |1⟩: 
            a1 = (control_bit == 0 ? 1.0 : -1.0) / sqrt(2.0) * phase_factor
            
            # Output states: copy bits change, main bits stay the same
            # For hadamard output |0⟩: copy targets unchanged from input
            copy_out_0 = copy(target_bits)
            push!(copy_out_0, 0)
            
            # For hadamard output |1⟩: copy targets unchanged from input (phase doesn't change basis)
            copy_out_1 = copy(target_bits)
            push!(copy_out_1, 1)
            
            v_exp = zeros(ComplexF64, 1 << (2k))
            
            # State with hadamard output |0⟩
            idx_0 = 0
            for l in 1:k
                idx_0 += main_bits[l] * (1 << (2*(l-1)))
                idx_0 += copy_out_0[l] * (1 << (2*(l-1) + 1))
            end
            v_exp[idx_0 + 1] = a0
            
            # State with hadamard output |1⟩
            idx_1 = 0
            for l in 1:k
                idx_1 += main_bits[l] * (1 << (2*(l-1)))
                idx_1 += copy_out_1[l] * (1 << (2*(l-1) + 1))
            end
            v_exp[idx_1 + 1] = a1
            
            @test isapprox(v_out, v_exp; atol=1e-10)
        end
    end
end
