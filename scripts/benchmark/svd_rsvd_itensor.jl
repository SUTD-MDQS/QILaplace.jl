# scripts/benchmark/svd_rsvd_itensor.jl
#
# Fig A runner — benchmark `svd(A, Linds...)` vs `QILaplace.RSVD.rsvd(A, Linds...)`
# on a random ITensor of size 2^n with the middle-site bipartition.
#
# Writes scripts/benchmark/results/svd_rsvd_itensor.jld2.
#
# Each method has its own time budget and cache / re-run switches (same layout
# idea as `qft_vs_fftw.jl`). SVD hitting its limit does not stop RSVD from
# continuing to larger n, and vice versa.

using ITensors
using LinearAlgebra
using Random

include(joinpath(@__DIR__, "common.jl"))

using QILaplace.RSVD: rsvd

const ART_PATH = joinpath(RESULTS_DIR, "svd_rsvd_itensor.jld2")

# ------------------------------- sweep / params ----------------------------
const N_RANGE      = 4:2:30
const RSVD_K       = 100
const RSVD_P       = 5
const RSVD_Q       = 2

const TIME_TO_STOP_SVD  = 10.0
const TIME_TO_STOP_RSVD = 10.0

const REBENCHMARK_SVD  = false
const REBENCHMARK_RSVD = false
const CACHE_ONLY_SVD   = false
const CACHE_ONLY_RSVD  = false

function load_state()
    data = load_results(ART_PATH)
    svd_s  = haskey(data, "svd")  ? from_dict(data["svd"])  : Series()
    rsvd_s = haskey(data, "rsvd") ? from_dict(data["rsvd"]) : Series()
    return svd_s, rsvd_s
end

function save_state(svd_s::Series, rsvd_s::Series)
    meta = machine_meta(; params = Dict(
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

function make_tensor(n::Int; seed::Int = 1234)
    rng = MersenneTwister(seed + n)
    sites = [Index(2, "site-$i") for i in 1:n]
    ITensors.disable_warn_order()
    A = randomITensor(sites...)
    return A, sites
end

function middle_inds(sites)
    n = length(sites)
    return sites[1:cld(n, 2)]
end

function run()
    banner("Fig A: SVD vs RSVD on a random ITensor")
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
        println("\n--- n = $n (matrix 2^$(cld(n, 2)) × 2^$(n - cld(n, 2))) ---")

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

        A, sites = make_tensor(n)
        Linds = middle_inds(sites)

        if exec_rsvd
            println("  RSVD: benchmarking (k=$RSVD_K, p=$RSVD_P, q=$RSVD_Q) ...")
            _ = rsvd(A, Linds...; k = RSVD_K, p = RSVD_P, q = RSVD_Q, verbose = false)
            trial = @benchmark rsvd($A, $(Linds)...; k = $RSVD_K, p = $RSVD_P, q = $RSVD_Q, verbose = false) samples = 5 evals = 1
            stats = stats_from_trial(trial)
            update!(rsvd_s, n, stats)
            @printf("  RSVD: time=%.4f s  mem=%.2f MB  allocs=%d  gc=%.4f s\n",
                    stats.time, stats.mem, stats.allocs, stats.gctime)
            save_state(svd_s, rsvd_s)
            if stats.time > TIME_TO_STOP_RSVD
                skip_rsvd = true
                println("  RSVD: exceeded TIME_TO_STOP_RSVD=$(TIME_TO_STOP_RSVD)s; larger n will not run for RSVD.")
            end
        end

        if exec_svd
            println("  SVD: benchmarking full svd ...")
            _ = svd(A, Linds...)
            trial = @benchmark svd($A, $(Linds)...) samples = 5 evals = 1
            stats = stats_from_trial(trial)
            update!(svd_s, n, stats)
            @printf("  SVD: time=%.4f s  mem=%.2f MB  allocs=%d  gc=%.4f s\n",
                    stats.time, stats.mem, stats.allocs, stats.gctime)
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
