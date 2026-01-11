using ITensors, Test
using QILaplace
using QILaplace.SignalConverters: _array_to_tensor, _tensor_to_mps_svd, _tensor_to_mps_rsvd, _tensor_to_mps

import LinearAlgebra.norm as la_norm
# using DataStructures: OrderedDict

# # Legacy function for comparison purposes
# function build_2n_mps_from_n(n::Int, mps::Vector{ITensor}, c_indices::OrderedDict{String, ITensors.Index}; threshold=1e-10)
#     sites = collect(values(c_indices))
#     paired_mps = Vector{ITensor}(undef, n)

#     Inds = inds(mps[1])
#     b_right = Inds[2]
#     paired_mps[1] = ITensor(sites[1], sites[2], b_right)
#     for br in 1:dim(b_right)
#         paired_mps[1][1, 1, br] = mps[1][1, br]
#         paired_mps[1][2, 2, br] = mps[1][2, br]
#     end

#     for i in 2:n-1
#         Inds = inds(mps[i])
#         b_left =  Inds[1]
#         b_right = Inds[3]
#         paired_mps[i] = ITensor(b_left, sites[2i-1], sites[2i], b_right)
#         for bl in 1:dim(b_left)
#             for br in 1:dim(b_right)
#                 paired_mps[i][bl, 1, 1, br] = mps[i][bl, 1, br]
#                 paired_mps[i][bl, 2, 2, br] = mps[i][bl, 2, br]
#             end 
#         end
#     end

#     Inds = inds(mps[n])
#     b_left =  Inds[1]
#     paired_mps[n] = ITensor(b_left, sites[2n-1], sites[2n])
#     for bl in 1:dim(b_left)
#         paired_mps[n][bl, 1, 1] = mps[n][bl, 1]
#         paired_mps[n][bl, 2, 2] = mps[n][bl, 2]
#     end

#     pairedMPS = Vector{ITensor}(undef, 2n)
#     U, S, V = svd(paired_mps[1], sites[1]; cutoff=threshold)
#     pairedMPS[1] = U
#     pairedMPS[2] =  S * V

#     for i in 2:n
#         Inds = inds(paired_mps[i])
#         U, S, V = svd(paired_mps[i], (Inds[1], Inds[2]); cutoff=threshold)
#         pairedMPS[2*i-1] = U
#         pairedMPS[2*i] =  S * V
#     end

#     return pairedMPS
# end


@testset "tensor_to_mps: RSVD vs SVD correctness" begin
    n = 6
    x = randn(2^n)
    ψ, _ = QILaplace.SignalConverters._array_to_tensor(x)

    mps_svd = QILaplace.SignalConverters._tensor_to_mps_svd(ψ; cutoff=1e-12, maxdim=128)
    mps_rsvd = QILaplace.SignalConverters._tensor_to_mps_rsvd(ψ; cutoff=1e-12, maxdim=128, k=64, p=5, q=2)

    # Reconstruct full tensors
    T_svd = contract_chain(mps_svd.data)
    T_rsvd = contract_chain(mps_rsvd.data)

    arr_orig = Array(ψ, inds(ψ)...)
    arr_svd = Array(T_svd, inds(ψ)...)
    arr_rsvd = Array(T_rsvd, inds(ψ)...)

    # Check approximate equality
    @test isapprox(arr_orig, arr_svd; atol=1e-9, rtol=0)
    @test isapprox(arr_orig, arr_rsvd; atol=1e-8, rtol=0)

    # Bond counts
    @test length(mps_svd.bonds) == n - 1
    @test length(mps_rsvd.bonds) == n - 1
end

@testset "signal_ztmps: 2-qubit copied-register identity" begin
    x = ComplexF64[
        1.0 + 0.0im,
        2.0 - 1.0im,
        -0.5 + 0.3im,
        0.25 - 0.75im,
    ]

    ψ_signal, α = QILaplace.signal_mps(x; cutoff=0.0, maxdim=typemax(Int))
    @test isapprox(la_norm(dense_from_signal_mps(ψ_signal)), 1.0; rtol=1e-12, atol=1e-12)

    zψ, α2 = QILaplace.signal_ztmps(x; cutoff=0.0, maxdim=typemax(Int))
    @test isapprox(α2, α; rtol=1e-12, atol=1e-12)

    ψ2n = QILaplace.Mps._as_signal_2n(zψ)
    A_sig = dense_from_signal_mps(ψ_signal)          # size (2,2)
    A_2n  = dense_from_signal_mps(ψ2n)               # size (2,2,2,2) in (m1,c1,m2,c2) order

    # Check: A_2n[m1,c1,m2,c2] = A_sig[m1,m2] if m1==c1 and m2==c2, else 0.
    for m1 in 1:2, c1 in 1:2, m2 in 1:2, c2 in 1:2
        if (m1 == c1) && (m2 == c2)
            @test isapprox(A_2n[m1, c1, m2, c2], A_sig[m1, m2]; rtol=1e-12, atol=1e-12)
        else
            @test isapprox(A_2n[m1, c1, m2, c2], 0.0 + 0.0im; rtol=0, atol=1e-12)
        end
    end
end

@testset "signal_ztmps vs legacy build_2n_mps_from_n (same 2n basis tensor)" begin
    x = ComplexF64[
        1.0 + 0.0im,
        2.0 - 1.0im,
        -0.5 + 0.3im,
        0.25 - 0.75im,
    ]

    ψ_signal, _ = QILaplace.signal_mps(x; cutoff=0.0, maxdim=typemax(Int))
    zψ, _ = QILaplace.signal_ztmps(x; cutoff=0.0, maxdim=typemax(Int))
    ψ2n_new = QILaplace.Mps._as_signal_2n(zψ)

    sites2n = ψ2n_new.sites
    c_indices = OrderedDict{String,Index}()
    for (i, s) in enumerate(sites2n)
        c_indices[string(i)] = s
    end

    pairedMPS_old = build_2n_mps_from_n(2, ψ_signal.data, c_indices; threshold=0.0)
    ψ2n_old = QILaplace.Mps.SignalMPS(pairedMPS_old, sites2n, 
        [commonind(pairedMPS_old[i], pairedMPS_old[i+1]) for i in 1:3])

    A_new = dense_from_signal_mps(ψ2n_new)
    A_old = dense_from_signal_mps(ψ2n_old)

    @test isapprox(A_old, A_new; rtol=1e-12, atol=1e-12)
end

@testset "Mps core invariants + adapters (round-trip)" begin
    x = ComplexF64[
        1.0 + 0.0im,
        2.0 + 0.0im,
        3.0 + 0.0im,
        4.0 + 0.0im,
    ]

    zψ, _ = QILaplace.signal_ztmps(x; cutoff=0.0, maxdim=typemax(Int))
    ψ2n = QILaplace.Mps._as_signal_2n(zψ)
    zψ_back = QILaplace.Mps._writeback_signal_2n(ψ2n)

    A0 = dense_from_signal_mps(QILaplace.Mps._as_signal_2n(zψ))
    A1 = dense_from_signal_mps(QILaplace.Mps._as_signal_2n(zψ_back))
    @test isapprox(A0, A1; rtol=1e-12, atol=1e-12)
end

@testset "norm correctness (SignalMPS vs dense contraction; zTMPS vs its 2n view)" begin
    x = randn(ComplexF64, 2^4)  # n=4
    ψ, _ = QILaplace.signal_mps(x; cutoff=0.0, maxdim=typemax(Int))
    T = contract_chain(ψ.data)
    ref = sqrt(abs(scalar(dag(T) * T)))
    @test isapprox(QILaplace.Mps.norm(ψ), ref; rtol=1e-10, atol=1e-10)

    zψ, _ = QILaplace.signal_ztmps(x; cutoff=0.0, maxdim=typemax(Int))
    @test isapprox(QILaplace.Mps.norm(zψ), QILaplace.Mps.norm(QILaplace.Mps._as_signal_2n(zψ)); rtol=1e-10, atol=1e-10)
end

@testset "compress! reduces bond dims (small sanity)" begin
    x = randn(ComplexF64, 2^6)  # n=6
    ψ, _ = QILaplace.signal_mps(x; cutoff=0.0, maxdim=typemax(Int))

    QILaplace.Mps.compress!(ψ; maxdim=4, tol=1e-12, sweeps=1)
    @test maximum(dim.(ψ.bonds)) ≤ 4

    zψ, _ = QILaplace.signal_ztmps(x; cutoff=0.0, maxdim=typemax(Int))
    QILaplace.Mps.compress!(zψ; maxdim=4, tol=1e-12, sweeps=1)
    @test maximum(dim.(zψ.bonds_main)) ≤ 4
    @test maximum(dim.(zψ.bonds_copy)) ≤ 4
end

@testset "update_bond! on zTMPS keeps PairCore.c consistent" begin
    x = randn(ComplexF64, 2^4)
    zψ, _ = QILaplace.signal_ztmps(x; cutoff=0.0, maxdim=typemax(Int))

    c_old = zψ.bonds_copy[2]
    c_new = Index(dim(c_old), "bond-copy-2-updated")

    QILaplace.Mps.update_bond!(zψ, c_old, c_new)

    @test zψ.bonds_copy[2] == c_new
    @test zψ.data[2].c == c_new
end
