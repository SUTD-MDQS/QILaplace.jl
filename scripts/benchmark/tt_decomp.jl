# scripts/benchmark/tt_decomp.jl
#
# Fig B runner — benchmark the full TT (MPS) decomposition of a random
# length-2^n signal via `signal_mps(x; method=:svd, ...)` versus
# `signal_mps(x; method=:rsvd, ...)`.
#
# Writes scripts/benchmark/results/tt_decomp.jld2.
#
# Controls (same spirit as `qft_vs_fftw.jl`): each method has its own time
# budget and its own cache / re-run switches. Exceeding the budget for SVD does
# not stop RSVD from continuing to larger n, and vice versa.

using ITensors
using LinearAlgebra
using Random

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "tt_decomp.jld2")

# ------------------------------- sweep / params ----------------------------
const N_RANGE      = 6:2:30
const CUTOFF       = 1e-12
const MAXDIM       = 1024
const RSVD_K       = 50
const RSVD_P       = 5
const RSVD_Q       = 2

# Mean @benchmark trial time (s) — per method; one method tripping its limit
# does not cancel the other series.
const TIME_TO_STOP_SVD  = 10.0
const TIME_TO_STOP_RSVD = 10.0

"""If true, start with an empty in-memory series so every n ∈ N_RANGE is re-benchmarked (artifact still merged on save)."""
const REBENCHMARK_SVD  = false
const REBENCHMARK_RSVD = false

"""If true, never run new work for that method — only report cached points or gaps."""
const CACHE_ONLY_SVD  = false
const CACHE_ONLY_RSVD = false

function load_state()
    data   = load_results(ART_PATH)
    svd_s  = haskey(data, "svd")  ? from_dict(data["svd"])  : Series()
    rsvd_s = haskey(data, "rsvd") ? from_dict(data["rsvd"]) : Series()
    return svd_s, rsvd_s
end

function save_state(svd_s::Series, rsvd_s::Series)
    meta = machine_meta(; params = Dict(
        "cutoff"  => CUTOFF,
        "maxdim"  => MAXDIM,
        "k"       => RSVD_K,
        "p"       => RSVD_P,
        "q"       => RSVD_Q,
        "n_range" => collect(N_RANGE),
        "time_to_stop_svd"  => TIME_TO_STOP_SVD,
        "time_to_stop_rsvd" => TIME_TO_STOP_RSVD,
        "rebenchmark_svd"    => REBENCHMARK_SVD,
        "rebenchmark_rsvd"   => REBENCHMARK_RSVD,
        "cache_only_svd"     => CACHE_ONLY_SVD,
        "cache_only_rsvd"    => CACHE_ONLY_RSVD,
    ))
    save_results(ART_PATH; meta = meta, series = Dict("svd" => svd_s, "rsvd" => rsvd_s))
    return nothing
end

function make_signal(n::Int; seed::Int = 1234)
    rng = MersenneTwister(seed + n)
    return randn(rng, 2^n)
end

function run()
    banner("Fig B: TT decomposition via signal_mps(:svd) vs signal_mps(:rsvd)")
    svd_s, rsvd_s = load_state()
    if REBENCHMARK_SVD
        svd_s = Series()
        println("REBENCHMARK_SVD=true — SVD series cleared in memory (will overwrite on save).")
    end
    if REBENCHMARK_RSVD
        rsvd_s = Series()
        println("REBENCHMARK_RSVD=true — RSVD series cleared in memory (will overwrite on save).")
    end

    skip_svd  = false
    skip_rsvd = false

    for n in N_RANGE
        println("\n--- n = $n (N = 2^$n = $(2^n)) ---")

        exec_rsvd = !CACHE_ONLY_RSVD && !skip_rsvd && (!has_point(rsvd_s, n) || REBENCHMARK_RSVD)
        exec_svd  = !CACHE_ONLY_SVD  && !skip_svd  && (!has_point(svd_s, n)  || REBENCHMARK_SVD)

        if CACHE_ONLY_RSVD
            println(has_point(rsvd_s, n) ?
                    "  RSVD: using cached point (CACHE_ONLY_RSVD=true)" :
                    "  RSVD: no cached point — CACHE_ONLY_RSVD=true; leaving n vacant")
        elseif skip_rsvd
            println("  RSVD: past time budget ($(TIME_TO_STOP_RSVD)s mean on an earlier n); not scheduling")
        elseif has_point(rsvd_s, n) && !REBENCHMARK_RSVD
            println("  RSVD: using cached point")
        end

        if CACHE_ONLY_SVD
            println(has_point(svd_s, n) ?
                    "  SVD: using cached point (CACHE_ONLY_SVD=true)" :
                    "  SVD: no cached point — CACHE_ONLY_SVD=true; leaving n vacant")
        elseif skip_svd
            println("  SVD: past time budget ($(TIME_TO_STOP_SVD)s mean on an earlier n); not scheduling")
        elseif has_point(svd_s, n) && !REBENCHMARK_SVD
            println("  SVD: using cached point")
        end

        if !exec_rsvd && !exec_svd
            continue
        end

        x = make_signal(n)

        if exec_rsvd
            println("  RSVD: benchmarking signal_mps(:rsvd) ...")
            ψ_warm, _ = signal_mps(x; method = :rsvd, cutoff = CUTOFF, maxdim = MAXDIM,
                                   k = RSVD_K, p = RSVD_P, q = RSVD_Q)
            mb = maxbond_mps(ψ_warm)
            trial = @benchmark signal_mps($x; method = :rsvd, cutoff = $CUTOFF, maxdim = $MAXDIM,
                                          k = $RSVD_K, p = $RSVD_P, q = $RSVD_Q) samples = 5 evals = 1
            stats = stats_from_trial(trial)
            update!(rsvd_s, n, stats; maxbond = mb)
            @printf("  RSVD: time=%.4f s  mem=%.2f MB  allocs=%d  maxbond=%d\n",
                    stats.time, stats.mem, stats.allocs, mb)
            save_state(svd_s, rsvd_s)
            if stats.time > TIME_TO_STOP_RSVD
                skip_rsvd = true
                println("  RSVD: exceeded TIME_TO_STOP_RSVD=$(TIME_TO_STOP_RSVD)s; larger n will not run for RSVD.")
            end
        end

        if exec_svd
            println("  SVD: benchmarking signal_mps(:svd) ...")
            ψ_warm, _ = signal_mps(x; method = :svd, cutoff = CUTOFF, maxdim = MAXDIM)
            mb = maxbond_mps(ψ_warm)
            trial = @benchmark signal_mps($x; method = :svd, cutoff = $CUTOFF, maxdim = $MAXDIM) samples = 5 evals = 1
            stats = stats_from_trial(trial)
            update!(svd_s, n, stats; maxbond = mb)
            @printf("  SVD: time=%.4f s  mem=%.2f MB  allocs=%d  maxbond=%d\n",
                    stats.time, stats.mem, stats.allocs, mb)
            save_state(svd_s, rsvd_s)
            if stats.time > TIME_TO_STOP_SVD
                skip_svd = true
                println("  SVD: exceeded TIME_TO_STOP_SVD=$(TIME_TO_STOP_SVD)s; larger n will not run for SVD.")
            end
        end
    end

    save_state(svd_s, rsvd_s)
    println("\nArtifact: ", ART_PATH)
    return nothing
end

run()
