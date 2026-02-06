# Tests for the full Z-Transform (ZT) transformer

ITensors.disable_warn_order()

using QILaplace.Mps: _as_signal_2n

# ----------------------------
# Analytical Reference
# ----------------------------

"""
    analytical_zt(x; ωr=2π, ωi=2π, Δt=1.0, normalize=true)

Compute the finite z-transform samples on the (k,ℓ) grid:
    s_{k,ℓ} = (ωr*k + im*ωi*ℓ)/N,   k,ℓ = 0..N-1
    χ[k+1,ℓ+1] = (normalize ? 1/N : 1) * Δt * ∑_{j=0}^{N-1} x[j+1] * exp(-s_{k,ℓ} * j * Δt)

Returns an N×N ComplexF64 matrix with canonical ordering (no bit reversal, no interleaving).
"""
function analytical_zt(
    x::AbstractVector; ωr::Real=2π, ωi::Real=2π, Δt::Real=1.0, normalize::Bool=true
)
    N = length(x)
    n = round(Int, log2(N))
    (1 << n) == N || throw(ArgumentError("length(x) must be a power of 2; got N=$N"))

    Z = Matrix{ComplexF64}(undef, N, N)
    pref = (normalize ? (1 / N) : 1.0) * Δt

    @inbounds for k in 0:(N - 1), ℓ in 0:(N - 1)
        s = (ωr * k + im * ωi * ℓ) / N
        acc = 0.0 + 0.0im
        for j in 0:(N - 1)
            acc += x[j + 1] * exp(-s * j * Δt)
        end
        Z[k + 1, ℓ + 1] = pref * acc
    end
    return Z
end

# ----------------------------
# Helper Functions
# ----------------------------

"""
    extract_zt_output(ψ_out::zTMPS) -> Matrix{ComplexF64}

Extract the Z-transform output matrix from a zTMPS for any n.
Returns a 2^n x 2^n matrix where rows correspond to main register (k)
and columns correspond to copy register (ℓ).
"""
function extract_zt_output(ψ_out::zTMPS)
    n = length(ψ_out.sites_main)
    ψ2n = _as_signal_2n(ψ_out)
    T = prod(ψ2n.data)

    # This maps to the (k, ℓ) grid where k is formed by main bits and ℓ by copy bits.
    ordered_inds = [ψ_out.sites_main; ψ_out.sites_copy]
    arr = Array(T, ordered_inds...)

    return reshape(ComplexF64.(arr), 2^n, 2^n) # implicitly assumes LSB first ordering
end

# ----------------------------
# Tests
# ----------------------------

@testset "zt_transformer.jl: Full zT on the basis states" begin
    ns = [1, 2, 3, 4]
    ωrs = [0.0, 0.75, 1.0, 2.0, 5.0]
    for (n, ωr) in Iterators.product(ns, ωrs)
        @testset "build_zt_mpo vs analytical_zt (n=$n)" begin
            N = 2^n
            ωi = 2π
            Δt = 1.0

            function bitreverse_int(x::Integer, n::Int)
                y = 0
                for i in 1:n
                    y <<= 1
                    y |= (x & 1)
                    x >>= 1
                end
                return y
            end

            for j in 0:(N - 1)
                x = zeros(Float64, N)
                x[j + 1] = 1.0

                # Build reference Z matrix
                Z_ref = analytical_zt(x; ωr=ωr, ωi=ωi, Δt=Δt, normalize=true)

                # Build zTMPS input and apply ZT MPO
                ψ_in, _ = signal_ztmps(x)
                mpo = build_zt_mpo(ψ_in, ωr)
                ψ_out = apply(mpo, ψ_in)

                # Extract output Z matrix
                Z_mpo = extract_zt_output(ψ_out)

                # Evaluate error
                err = LinearAlgebra.norm(Z_mpo - Z_ref)
                @test err ≤ 1e-7
            end
        end
    end
end
