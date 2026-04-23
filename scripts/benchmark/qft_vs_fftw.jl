# scripts/benchmark/qft_vs_fftw.jl
#
# Fig D runner — runtime comparison for:
# - QFT only      : apply(W, ψ) with pre-built (R)SVD MPS (solid in plot)
# - (R)SVD + QFT    : signal_mps(method=:rsvd) + apply(W, ψ) (dashed in plot)
# - FFTW baseline : a single random-data curve (black diamonds)
#
# Writes scripts/benchmark/results/qft_vs_fftw.jld2.

using FFTW, ITensors
using LinearAlgebra, Random, Printf

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "qft_vs_fftw.jld2")

# ------------------------------- sweep / params ----------------------------
const N_RANGE      = 4:2:28
const SIGNAL_KINDS = [:random, :sin, :sine20, :sin_cusp]
const CUTOFF       = 1e-12
const TIME_TO_STOP = 1.0
const RSVD_K       = 15
const MAXDIM       = 512
const RSVD_P       = 5
const RSVD_Q       = 2
# Toggle this single line to switch compression backend:
# const COMP_METHOD  = :svd
const COMP_METHOD  = :rsvd
const COMP_LABEL   = COMP_METHOD == :svd ? "SVD" : "RSVD"

# For an n-qubit signal reshaped as (2,2,...,2), the largest possible
# Schmidt rank across the middle bipartition is 2^(floor(n/2)).
max_middle_rank(n::Int) = 2 ^ fld(n, 2)

# Structured signals keep a fixed reduced rank; random signals use the
# maximal middle rank so we do not artificially cap bond growth.
rsvd_k_for_signal(kind::Symbol, n::Int) = kind === :random ? max_middle_rank(n) : RSVD_K

comp_kwargs(kind::Symbol, n::Int) = COMP_METHOD === :rsvd ?
    (; maxdim = rsvd_k_for_signal(kind, n), method = :rsvd, k = rsvd_k_for_signal(kind, n), p = RSVD_P, q = RSVD_Q) :
    (; method = :svd)

function load_state()
    data = load_results(ART_PATH)
    series = Dict{Symbol, Series}()
    for kind in SIGNAL_KINDS
        series[Symbol("qft_only_$(kind)")] = haskey(data, "qft_only_$(kind)") ? from_dict(data["qft_only_$(kind)"]) : Series()
        series[Symbol("rsvd_qft_$(kind)")] = haskey(data, "rsvd_qft_$(kind)") ? from_dict(data["rsvd_qft_$(kind)"]) : Series()
    end
    series[:fftw_random] = haskey(data, "fftw_random") ? from_dict(data["fftw_random"]) : Series()
    return series
end

function save_state(series::Dict{Symbol, Series})
    meta = machine_meta(; params = Dict(
        "comp_method" => String(COMP_METHOD),
        "cutoff"  => CUTOFF,
        "maxdim"  => MAXDIM,
        "rsvd_k_structured"  => RSVD_K,
        "rsvd_k_random"      => "2^(floor(n/2))",
        "rsvd_p"  => RSVD_P,
        "rsvd_q"  => RSVD_Q,
        "n_range" => collect(N_RANGE),
        "signal_kinds" => String.(SIGNAL_KINDS),
    ))
    out = Dict{String, Any}()
    for (k, s) in series
        out[String(k)] = s
    end
    save_results(ART_PATH; meta = meta, series = out)
    return nothing
end

qft_fftw(x::AbstractVector) = begin
    N = length(x)
    xhat = x ./ norm(x)
    bfft(xhat) ./ sqrt(N)
end

function run()
    banner("Fig D: QFT only vs $(COMP_LABEL)+QFT + FFTW(random)")
    series = load_state()
    skip = Dict{Symbol, Bool}(k => false for k in keys(series))

    for n in N_RANGE
        println("\n--- n = $n (N = $(2^n)) ---")

        pending_qft = any(kind -> begin
            key_qft_only = Symbol("qft_only_$(kind)")
            key_comp_qft = Symbol("rsvd_qft_$(kind)")
            needs_point = !has_point(series[key_qft_only], n) || !has_point(series[key_comp_qft], n)
            can_run = !(skip[key_qft_only] || skip[key_comp_qft])
            return needs_point && can_run
        end, SIGNAL_KINDS)

        pending_fftw = !has_point(series[:fftw_random], n) && !skip[:fftw_random]

        if !(pending_qft || pending_fftw)
            println("nothing pending at this n; skipping")
            continue
        end

        for kind in SIGNAL_KINDS
            key_qft_only = Symbol("qft_only_$(kind)")
            key_rsvd_qft = Symbol("rsvd_qft_$(kind)")
            need_qft_kind = (!has_point(series[key_qft_only], n) || !has_point(series[key_rsvd_qft], n)) &&
                            !(skip[key_qft_only] || skip[key_rsvd_qft])
            if !need_qft_kind
                println("  $kind already benchmarked; skipping")
                continue
            end
            
            # Run if we don't have the QFT only or (R)SVD+QFT benchmark data for this n and kind.
            x = make_signal(kind, n)
            ksig = COMP_METHOD === :rsvd ? rsvd_k_for_signal(kind, n) : 0

            # 1) Compression-only timing on this exact (n, kind) signal.
            ψ_apply, _ = signal_mps(x; cutoff = CUTOFF, comp_kwargs(kind, n)...)
            d_qft_only = maxbond_mps(ψ_apply)
            trial_comp = @benchmark signal_mps($x; cutoff = $CUTOFF, comp_kwargs($kind, $n)...) samples = 5 evals = 1
            stats_comp = stats_from_trial(trial_comp)

            W = build_qft_mpo(ψ_apply; cutoff = CUTOFF, maxdim = MAXDIM)

            # 2) Apply-only timing using an MPS built from the same data instance.
            _ = apply(W, ψ_apply) # dry run
            trial_apply = @benchmark apply($W, $ψ_apply) samples = 5 evals = 1
            stats_apply = stats_from_trial(trial_apply)

            # QFT only = apply-only.
            update!(series[key_qft_only], n, stats_apply; maxbond = d_qft_only)

            # (R)SVD+QFT = compression + apply.
            stats_total = RunStats(
                stats_comp.time + stats_apply.time,
                stats_comp.gctime + stats_apply.gctime,
                stats_comp.mem + stats_apply.mem,
                stats_comp.allocs + stats_apply.allocs,
            )
            update!(series[key_rsvd_qft], n, stats_total; maxbond = d_qft_only)

            @printf("  QFT only  [%s]: t_apply=%.4f s  mem=%.2f MB  Dmax=%d\n",
                    string(kind), stats_apply.time, stats_apply.mem, d_qft_only)
            if COMP_METHOD === :rsvd
                @printf("  %s+QFT  [%s]: t_comp=%.4f s  t_apply=%.4f s  total=%.4f s  Dmax=%d  k=%d\n",
                        COMP_LABEL, string(kind), stats_comp.time, stats_apply.time, stats_total.time, d_qft_only, ksig)
            else
                @printf("  %s+QFT  [%s]: t_comp=%.4f s  t_apply=%.4f s  total=%.4f s  Dmax=%d\n",
                        COMP_LABEL, string(kind), stats_comp.time, stats_apply.time, stats_total.time, d_qft_only)
            end
            save_state(series)
            if stats_total.time > TIME_TO_STOP
                skip[key_qft_only] = true
                skip[key_rsvd_qft] = true
                println("  qft_only[$kind] past time_to_stop; skipping larger n.")
            end
        end

        if has_point(series[:fftw_random], n) && !skip[:fftw_random]
            continue
        end
        
        # Run if we don't have a FFTW benchmark point for this n.
        xrand = make_signal(:random, n)
        _ = qft_fftw(xrand)
        trial = @benchmark qft_fftw($xrand) samples = 5 evals = 1
        stats = stats_from_trial(trial)
        update!(series[:fftw_random], n, stats)
        @printf("  FFTW      [random]: time=%.4f s  mem=%.2f MB\n", stats.time, stats.mem)
        save_state(series)
        if stats.time > TIME_TO_STOP
            skip[:fftw_random] = true
            println("  fftw_random past time_to_stop; skipping larger n.")
        end
    end

    save_state(series)
    println("\nArtifact: ", ART_PATH)
    return nothing
end

run()
