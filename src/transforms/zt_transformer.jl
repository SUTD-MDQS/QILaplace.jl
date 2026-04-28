#src/transforms/zt_transformer.jl
# This module builds the discrete Laplace transform aka z-transform that acts on a given MPS to apply the damping transform circuit. 
# It used the dt_gates from DTGates and zt_gates from ZTGates to build the full zT MPO as a PairedSiteMPO that acts on ZTMPS representations of signals. 

# The zT circuit is the composition of the damping transform (DT) and the paired 2n-site QFT
# used in ZTMPS (`control_Hphase_ztmps_mpo`). We build `W_dt` and `W_qft` separately, then
# fuse once (`apply(W_dt, W_qft)`) and run a final compress sweep

module ZTTransformer

using ITensors
using ..Mpo: PairedSiteMPO
using ..Mps: ZTMPS
using ..DTGates: control_damping_mpo, control_damping_copy_mpo
using ..ZTGates: control_Hphase_ztmps_mpo
using ..DTTransform: build_dt_mpo, zip_to_combine_mpos, zip_to_compress_mpo
using ..ApplyMPO: apply

export build_zt_mpo

"""
    build_zt_mpo(n, ωr, sites_main, sites_copy; cutoff=1e-14, maxdim=1000) -> PairedSiteMPO
    build_zt_mpo(ψ::ZTMPS, ωr; cutoff=1e-14, maxdim=1000)

Build the discrete Laplace transform (z-Transform) MPO for `n` qubits at
damping parameter `ωr`. The cutoff tells the error threshold of the compressed MPO.

Operationally, this builds the **DT** (`build_dt_mpo`) and the **paired ZTMPS QFT**
(`control_Hphase_ztmps_mpo` on all `2n` sites), then fuses them with `apply(W_dt, W_qft)` and
compresses — same end operator as the legacy interleaved construction, without per-stage DT×QFT growth.

When called with a `ZTMPS`, site indices are taken from the MPS itself.

# Examples
```julia
ψ_z = signal_ztmps(x)
W_zt = build_zt_mpo(ψ_z, ωr)
ψ_out = apply(W_zt, ψ_z)
```
"""
function build_zt_mpo(
    n::Int,
    ωr::Real,
    sites_main::Vector{I},
    sites_copy::Vector{I};
    cutoff=1e-14,
    maxdim=1000,
) where {I<:Index}
    n >= 1 || throw(ArgumentError("build_zt_mpo: n must be ≥ 1. Found n=$n"))
    length(sites_main) == n || throw(
        ArgumentError(
            "build_zt_mpo: sites_main must have $n elements. Got $(length(sites_main))"
        ),
    )
    length(sites_copy) == n || throw(
        ArgumentError(
            "build_zt_mpo: sites_copy must have $n elements. Got $(length(sites_copy))"
        ),
    )

    # Interleave sites: [main[1], copy[1], main[2], copy[2], ...]
    all_sites = Vector{I}(undef, 2n)
    for i in 1:n
        all_sites[2i-1] = sites_main[i]
        all_sites[2i] = sites_copy[i]
    end

    if n == 1
        W_dt = build_dt_mpo(1, ωr, sites_main, sites_copy; cutoff=cutoff, maxdim=maxdim)
        W_qft = control_Hphase_ztmps_mpo(1, all_sites)
        return apply(W_dt, W_qft)
    end

    # ---- DT MPO (damping only), same sites as before ----
    W_dt = build_dt_mpo(n, ωr, sites_main, sites_copy; cutoff=cutoff, maxdim=maxdim)

    # ---- Full paired 2n-site QFT MPO (zT QFT blocks only), then fuse with DT once ----
    mpo_qft = control_Hphase_ztmps_mpo(1, all_sites[1:2])
    oc_q = 0
    for k in 2:n
        if length(mpo_qft.sites_main) < k
            s_main = all_sites[2*k-1]
            s_copy = all_sites[2*k]
            b_main = Index(1, "bond-main-$(k-1)")
            b_copy = Index(1, "bond-copy-$k")
            mpo_qft.data[end] *= ITensor(1.0, b_main)
            push!(mpo_qft.bonds_main, b_main)
            T_main = delta(s_main, s_main') * ITensor(1.0, b_main) * ITensor(1.0, b_copy)
            T_copy = delta(s_copy, s_copy') * ITensor(1.0, b_copy)
            push!(mpo_qft.data, T_main)
            push!(mpo_qft.data, T_copy)
            push!(mpo_qft.sites_main, s_main)
            push!(mpo_qft.sites_copy, s_copy)
            push!(mpo_qft.bonds_copy, b_copy)
        end
        block_k = control_Hphase_ztmps_mpo(k, all_sites[1:2k])
        mpo_qft, oc_q, _ = zip_to_combine_mpos(mpo_qft, block_k, oc_q)
        mpo_qft, oc_q = zip_to_compress_mpo(mpo_qft, oc_q, "down"; cutoff=cutoff, maxdim=maxdim)
    end

    # Paired MPO composition must match the legacy interleaved zip-up; empirically this is `apply(W_dt, W_qft)`
    # (see `ApplyMPO` `T1 * T2` layout), not `apply(W_qft, W_dt)`.
    W_zt = apply(W_dt, mpo_qft)
    W_zt, _ = zip_to_compress_mpo(W_zt, 1, "down"; cutoff=cutoff, maxdim=maxdim)
    return W_zt
end

function build_zt_mpo(ψ::ZTMPS, ωr::Real; cutoff=1e-14, maxdim=1000)
    build_zt_mpo(
        length(ψ.sites_main), ωr, ψ.sites_main, ψ.sites_copy; cutoff=cutoff, maxdim=maxdim
    )
end

end # module ZTTransformer
