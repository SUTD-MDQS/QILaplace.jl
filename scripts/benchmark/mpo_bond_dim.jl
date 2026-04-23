# scripts/benchmark/mpo_bond_dim.jl
#
# Fig C runner — build the QFT, DT, and zT MPOs for each qubit count `n` in a
# sweep and record their bond dimension proxy from the compressed MPO bond
# indices (equivalent to retained rank after cutoff). This reproduces the
# bounded-bond-dimension story from the paper (left panel of the deprecated
# plot_benchmark_time.jl).
#
# Writes scripts/benchmark/results/mpo_bond_dim.jld2.

using ITensors
using LinearAlgebra
using Random

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "mpo_bond_dim.jld2")

# ------------------------------- sweep / params ----------------------------
const N_RANGE = 2:1:30           # qft sites = n; dt/zt operate on 2n sites
const OMEGA_R = 2π
const CUTOFF  = 1e-15
const MAXDIM  = typemax(Int)
const BOND_METRIC = "maxbond_mpo_after_cutoff"
const TIME_TO_STOP = 30.0

function load_state()
    data = load_results(ART_PATH)
    if haskey(data, "meta")
        params = get(data["meta"], "params", Dict{String, Any}())
        old_cutoff = get(params, "cutoff", nothing)
        old_n_range = get(params, "n_range", nothing)
        old_metric = get(params, "bond_metric", nothing)
        if old_cutoff != CUTOFF || old_n_range != collect(N_RANGE) || old_metric != BOND_METRIC
            @warn "Existing artifact metadata does not match current legacy-aligned settings; starting fresh." old_cutoff old_n_range old_metric
            return Series(), Series(), Series()
        end
    end
    get_ser(k) = haskey(data, k) ? from_dict(data[k]) : Series()
    return get_ser("qft"), get_ser("dt"), get_ser("zt")
end

function save_state(qft_s::Series, dt_s::Series, zt_s::Series)
    meta = machine_meta(; params = Dict(
        "cutoff"  => CUTOFF,
        "maxdim"  => MAXDIM,
        "omega_r" => OMEGA_R,
        "n_range" => collect(N_RANGE),
        "bond_metric" => BOND_METRIC,
        "time_to_stop" => TIME_TO_STOP,
    ))
    save_results(ART_PATH; meta = meta, series = Dict(
        "qft" => qft_s, "dt" => dt_s, "zt" => zt_s,
    ))
    return nothing
end

function bench_qft_mpo(n::Int)
    x = make_signal(:sin, n)
    ψ, _ = signal_mps(x; cutoff = CUTOFF, maxdim = MAXDIM)
    t0 = time()
    W = build_qft_mpo(ψ; cutoff = CUTOFF, maxdim = MAXDIM)
    elapsed = time() - t0
    return maxbond_mpo(W), elapsed
end

function bench_dt_mpo(n::Int)
    x = make_signal(:sin, n)
    ψ_z, _ = signal_ztmps(x; cutoff = CUTOFF, maxdim = MAXDIM)
    t0 = time()
    W = build_dt_mpo(ψ_z, OMEGA_R; cutoff = CUTOFF, maxdim = MAXDIM)
    elapsed = time() - t0
    return maxbond_mpo(W), elapsed
end

function bench_zt_mpo(n::Int)
    x = make_signal(:sin, n)
    ψ_z, _ = signal_ztmps(x; cutoff = CUTOFF, maxdim = MAXDIM)
    t0 = time()
    W = build_zt_mpo(ψ_z, OMEGA_R; cutoff = CUTOFF, maxdim = MAXDIM)
    elapsed = time() - t0
    return maxbond_mpo(W), elapsed
end

function run()
    banner("Fig C: MPO bond dimension vs n for QFT / DT / zT")
    qft_s, dt_s, zt_s = load_state()
    skip_qft = false
    skip_dt = false
    skip_zt = false

    for n in N_RANGE
        println("\n--- n = $n ---")

        if !has_point(qft_s, n) && !skip_qft
            try
                mb, elapsed = bench_qft_mpo(n)
                update!(qft_s, n, RunStats(elapsed, 0.0, 0.0, 0); maxbond = mb)
                @printf("  QFT MPO: maxbond=%d  build=%.3f s\n", mb, elapsed)
                save_state(qft_s, dt_s, zt_s)
                if elapsed > TIME_TO_STOP
                    skip_qft = true
                    println("  QFT past time_to_stop; skipping larger n.")
                end
            catch err
                @warn "QFT MPO failed at n=$n" exception = err
                skip_qft = true
            end
        end

        if !has_point(dt_s, n) && !skip_dt
            try
                mb, elapsed = bench_dt_mpo(n)
                update!(dt_s, n, RunStats(elapsed, 0.0, 0.0, 0); maxbond = mb)
                @printf("  DT  MPO: maxbond=%d  build=%.3f s\n", mb, elapsed)
                save_state(qft_s, dt_s, zt_s)
                if elapsed > TIME_TO_STOP
                    skip_dt = true
                    println("  DT past time_to_stop; skipping larger n.")
                end
            catch err
                @warn "DT MPO failed at n=$n" exception = err
                skip_dt = true
            end
        end

        if !has_point(zt_s, n) && !skip_zt
            try
                mb, elapsed = bench_zt_mpo(n)
                update!(zt_s, n, RunStats(elapsed, 0.0, 0.0, 0); maxbond = mb)
                @printf("  zT  MPO: maxbond=%d  build=%.3f s\n", mb, elapsed)
                save_state(qft_s, dt_s, zt_s)
                if elapsed > TIME_TO_STOP
                    skip_zt = true
                    println("  zT past time_to_stop; skipping larger n.")
                end
            catch err
                @warn "zT MPO failed at n=$n" exception = err
                skip_zt = true
            end
        end
    end

    save_state(qft_s, dt_s, zt_s)
    println("\nArtifact: ", ART_PATH)
    return nothing
end

run()
