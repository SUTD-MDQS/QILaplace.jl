using QILaplace.RSVD: rsvd

const RSVD_REL_TOL = 1e-10

function low_rank_fixture(dim_left::Int, dim_right::Int, r_true::Int; seed::Int=42)
    Random.seed!(seed)
    i = Index(dim_left, "i")
    j = Index(dim_right, "j")
    bond_k = Index(r_true, "k")
    bond_l = Index(r_true, "l")
    U_true = random_itensor(i, bond_k)
    singulars = [exp(-0.5 * n) for n in 1:r_true]
    S_true = diag_itensor(singulars, bond_k, bond_l)
    V_true = random_itensor(bond_l, j)
    return (U_true * S_true * V_true, i, j, singulars)
end

function singular_values(S::ITensor)
    inds_S = inds(S)
    @assert length(inds_S) == 2
    left_s, right_s = inds_S
    Sd = Array(S, left_s, right_s)
    n = min(size(Sd, 1), size(Sd, 2))
    return [Sd[m, m] for m in 1:n]
end

@testset "rsvd.jl: correctness, parameter controls, error handling" begin
    @testset "Low-rank reconstruction" begin
        A, i, j, _ = low_rank_fixture(100, 100, 10)
        k = 15
        U, S, V = rsvd(A, i; k=k, p=5, q=2, cutoff=1e-10)

        A_approx = U * S * V
        rel_err = ITensors.norm(A - A_approx) / ITensors.norm(A)
        @test rel_err < RSVD_REL_TOL

        u_link = commonind(U, S)
        v_link = commonind(S, V)
        @test dim(u_link) <= k + 5

        U_iso = U * dag(prime(U, u_link))
        V_iso = dag(prime(V, v_link)) * V
        Id_u = diag_itensor(1.0, u_link, u_link')
        Id_v = diag_itensor(1.0, v_link', v_link)
        @test ITensors.norm(U_iso - Id_u) < 1e-10
        @test ITensors.norm(V_iso - Id_v) < 1e-10

        s_vals = singular_values(S)
        @test all(abs.(imag.(s_vals)) .< 1e-12)
        real_s = real.(s_vals)
        @test all(real_s .>= -1e-12)
        @test issorted(real_s; rev=true)
        Sdense = Array(S, inds(S)...)
        @test LinearAlgebra.norm(Sdense - Diagonal(real_s)) < 1e-10

        U_trunc, S_trunc, V_trunc = rsvd(A, i; k=5, p=5, q=2, maxdim=5)
        u_link_trunc = commonind(U_trunc, S_trunc)
        @test dim(u_link_trunc) <= 5

        U_smoke, S_smoke, V_smoke = rsvd(A, i; k=5, p=2, q=1)
        @test U_smoke isa ITensor && S_smoke isa ITensor && V_smoke isa ITensor
    end

    @testset "Parameter controls" begin
        A, i, _, _ = low_rank_fixture(40, 35, 6)

        cutoff = 5e-2
        U_lo, S_lo, _ = rsvd(A, i; k=10, p=4, cutoff=1e-12)
        U_cut, S_cut, _ = rsvd(A, i; k=10, p=4, cutoff=cutoff)
        kept_lo = dim(commonind(U_lo, S_lo))
        kept_hi = dim(commonind(U_cut, S_cut))
        @test kept_hi < kept_lo
        @test kept_hi ≥ 1

        mindim = 4
        U_min, S_min, _ = rsvd(A, i; k=10, p=2, cutoff=0.5, mindim=mindim, maxdim=10)
        @test dim(commonind(U_min, S_min)) ≥ mindim

        custom_tag = "Link,mybond"
        U_tag, S_tag, V_tag = rsvd(A, i; k=8, bondtag=custom_tag)
        link_left = commonind(U_tag, S_tag)
        link_right = commonind(S_tag, V_tag)
        @test occursin("mybond", string(tags(link_left)))
        @test occursin("mybond", string(tags(link_right)))

        U_seed1, S_seed1, V_seed1 = rsvd(A, i; k=8, p=3, random_seed=2024)
        U_seed2, S_seed2, V_seed2 = rsvd(A, i; k=8, p=3, random_seed=2024)
        A1 = U_seed1 * S_seed1 * V_seed1
        A2 = U_seed2 * S_seed2 * V_seed2
        @test ITensors.norm(A1 - A2) < 1e-12
    end

    @testset "Multiple left indices" begin
        i1 = Index(6, "i1")
        i2 = Index(5, "i2")
        j1 = Index(4, "j1")
        j2 = Index(3, "j2")
        bond = Index(7, "r")
        Random.seed!(123)
        U_true = random_itensor(i1, i2, bond)
        S_true = diag_itensor([0.9^n for n in 1:dim(bond)], bond, prime(bond))
        V_true = random_itensor(prime(bond), j1, j2)
        A = U_true * S_true * V_true

        U_multi, S_multi, V_multi = rsvd(A, i1, i2; k=7, p=2, q=1)
        A_hat = U_multi * S_multi * V_multi
        rel_err = ITensors.norm(A - A_hat) / ITensors.norm(A)
        @test rel_err < RSVD_REL_TOL
    end

    @testset "Error handling" begin
        i = Index(8, "i")
        j = Index(8, "j")
        A = random_itensor(i, j)
        bogus = Index(5, "bogus")

        @test_throws ErrorException rsvd(A; k=4)
        @test_throws ErrorException rsvd(A, bogus; k=4)
    end
end
