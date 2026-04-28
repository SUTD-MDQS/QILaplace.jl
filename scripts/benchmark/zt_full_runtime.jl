# scripts/benchmark/zt_full_runtime.jl
#
# Fig E runner — zT pipeline timing aligned with `qft_vs_fftw.jl`:
# - zT only      : `apply(W, ψ)` with W = build_zt_mpo(ψ, ωr) (solid in plot)
# - Encoding+zT  : `signal_ztmps` + `apply` on the same ψ instance (dashed)
#
# For each (n, signal kind): build ψ once for MPO construction, benchmark
# `signal_ztmps` (encoding), dry-run `apply`, then benchmark `apply(W, ψ)`.
# The zT MPO is built outside the timed encoding region (same separation idea
# as building `build_qft_mpo(ψ)` after the compression MPS in the QFT script).
#
# Signals: `:sin`, `:multi_sin_exp` (damped multi-sine),
# `:abs_cos_power_p8` (|cos|^{0.8}), `:random` (Gaussian noise, same as QFT).
#
# Writes scripts/benchmark/results/zt_full_runtime.jld2 (schema unchanged for
# `plot_zt_runtime.jl`: per-kind dicts with `t_input_mps`, `t_apply`, `maxbond`).

using BenchmarkTools
using ITensors
using LinearAlgebra
using Printf

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "zt_full_runtime.jld2")

# ------------------------------- sweep / params ----------------------------
const N_RANGE      = 2:1:30
const SIGNAL_KINDS = [:sin, :multi_sin_exp, :abs_cos_power_p8, :random]
const OMEGA_R      = 2π
const CUTOFF       = 1e-15
const MAXDIM       = 512
const TIME_TO_STOP = 100.0
const RSVD_K       = 15
const RSVD_P       = 5
const RSVD_Q       = 2
# Toggle compression backend for the inner `signal_mps` inside `signal_ztmps`:
# const COMP_METHOD  = :svd
const COMP_METHOD  = :rsvd
const COMP_LABEL   = COMP_METHOD == :svd ? "SVD" : "RSVD"

max_middle_rank(n::Int) = 2 ^ fld(n, 2)
rsvd_k_for_signal(kind::Symbol, n::Int) = kind === :random ? max_middle_rank(n) : RSVD_K

function zt_comp_kwargs(kind::Symbol, n::Int)
    if COMP_METHOD === :rsvd
        k = rsvd_k_for_signal(kind, n)
        return (; cutoff = CUTOFF, maxdim = k, method = :rsvd, k = k, p = RSVD_P, q = RSVD_Q)
    else
        return (; cutoff = CUTOFF, maxdim = MAXDIM, method = :svd)
    end
end

# ---------------------------------------------------------------------------
# On-disk layout (unchanged for plot_zt_runtime.jl)
# ---------------------------------------------------------------------------

mutable struct ZTSeries
    t_input_mps::Dict{Int, Float64}
    t_apply::Dict{Int, Float64}
    maxbond::Dict{Int, Int}
end

ZTSeries() = ZTSeries(Dict{Int, Float64}(), Dict{Int, Float64}(), Dict{Int, Int}())

has_point(s::ZTSeries, n::Int) = haskey(s.t_input_mps, n) && haskey(s.t_apply, n)

to_dict(s::ZTSeries) = Dict{String, Any}(
    "t_input_mps" => s.t_input_mps,
    "t_apply"     => s.t_apply,
    "maxbond"     => s.maxbond,
)

function from_dict_zt(d::AbstractDict)
    s = ZTSeries()
    haskey(d, "t_input_mps") && merge!(s.t_input_mps, d["t_input_mps"])
    haskey(d, "t_apply")     && merge!(s.t_apply,     d["t_apply"])
    haskey(d, "maxbond")     && merge!(s.maxbond,     d["maxbond"])
    return s
end

function meta_compatible(params::AbstractDict)::Bool
    get(params, "cutoff", nothing) == CUTOFF || return false
    get(params, "n_range", nothing) == collect(N_RANGE) || return false
    get(params, "signal_kinds", nothing) == String.(SIGNAL_KINDS) || return false
    get(params, "maxdim", nothing) == MAXDIM || return false
    get(params, "comp_method", nothing) == String(COMP_METHOD) || return false
    get(params, "rsvd_k_structured", nothing) == RSVD_K || return false
    return true
end

function load_state()
    data = load_results(ART_PATH)
    if haskey(data, "meta")
        params = get(data["meta"], "params", Dict{String, Any}())
        if !meta_compatible(params)
            @warn "Existing zt_full_runtime artifact params differ from this script; starting fresh." params
            return Dict(kind => ZTSeries() for kind in SIGNAL_KINDS)
        end
    end
    series = Dict{Symbol, ZTSeries}()
    for kind in SIGNAL_KINDS
        label = String(kind)
        series[kind] = haskey(data, label) ? from_dict_zt(data[label]) : ZTSeries()
    end
    return series
end

function save_state(series::Dict{Symbol, ZTSeries})
    meta = machine_meta(; params = Dict(
        "omega_r"      => OMEGA_R,
        "cutoff"       => CUTOFF,
        "maxdim"       => MAXDIM,
        "n_range"      => collect(N_RANGE),
        "signal_kinds" => String.(SIGNAL_KINDS),
        "comp_method"  => String(COMP_METHOD),
        "rsvd_k_structured" => RSVD_K,
        "rsvd_k_random"     => "2^(floor(n/2))",
        "rsvd_p"       => RSVD_P,
        "rsvd_q"       => RSVD_Q,
        "time_to_stop" => TIME_TO_STOP,
    ))
    out = Dict{String, Any}()
    for (k, s) in series
        out[String(k)] = to_dict(s)
    end
    save_results(ART_PATH; meta = meta, series = out)
    return nothing
end

function skip_dict()
    d = Dict{Symbol, Bool}()
    for kind in SIGNAL_KINDS
        d[Symbol("zt_only_$(kind)")]   = false
        d[Symbol("encode_zt_$(kind)")] = false
    end
    return d
end

const REBENCHMARK = false

function run()
    banner("Fig E: zT apply only vs $(COMP_LABEL)+apply (per signal kind)")
    series = load_state()
    if REBENCHMARK
        series = Dict(kind => ZTSeries() for kind in SIGNAL_KINDS)
        println("REBENCHMARK=true — cleared all series in memory.")
    end
    skip = skip_dict()

    for n in N_RANGE
        println("\n--- n = $n (output qubits m = $(2n)) ---")

        pending_zt = any(kind -> begin
            k1 = Symbol("zt_only_$(kind)")
            k2 = Symbol("encode_zt_$(kind)")
            needs = !has_point(series[kind], n)
            can_run = !(skip[k1] || skip[k2])
            return needs && can_run
        end, SIGNAL_KINDS)

        if !pending_zt
            println("nothing pending at this n; skipping")
            continue
        end

        for kind in SIGNAL_KINDS
            k1 = Symbol("zt_only_$(kind)")
            k2 = Symbol("encode_zt_$(kind)")
            need_kind = !has_point(series[kind], n) && !(skip[k1] || skip[k2])
            if !need_kind
                println("  $(string(kind)) already benchmarked; skipping")
                continue
            end

            x = make_signal(kind, n)
            ksig = COMP_METHOD === :rsvd ? rsvd_k_for_signal(kind, n) : 0

            # 1) Encoding MPS for site/bond structure + zT MPO (outside timed encoding).
            ψ_apply = signal_ztmps(x; zt_comp_kwargs(kind, n)...)
            d_zt = maxbond_mps(ψ_apply)
            W = build_zt_mpo(ψ_apply, OMEGA_R; cutoff = CUTOFF, maxdim = MAXDIM)

            # 2) Timed encoding (same kwargs as the warm construction).
            trial_comp = @benchmark signal_ztmps($x; zt_comp_kwargs($kind, $n)...) samples = 5 evals = 1
            stats_comp = stats_from_trial(trial_comp)

            # 3) Apply-only: same ψ and W as in the QFT script (apply does not mutate inputs).
            _ = apply(W, ψ_apply)
            trial_apply = @benchmark apply($W, $ψ_apply) samples = 5 evals = 1
            stats_apply = stats_from_trial(trial_apply)

            series[kind].t_input_mps[n] = stats_comp.time
            series[kind].t_apply[n]     = stats_apply.time
            series[kind].maxbond[n]     = d_zt

            stats_total = RunStats(
                stats_comp.time + stats_apply.time,
                stats_comp.gctime + stats_apply.gctime,
                stats_comp.mem + stats_apply.mem,
                stats_comp.allocs + stats_apply.allocs,
            )

            @printf("  zT only   [%s]: t_apply=%.4f s  mem=%.2f MB  Dmax=%d\n",
                    string(kind), stats_apply.time, stats_apply.mem, d_zt)
            if COMP_METHOD === :rsvd
                @printf("  %s+zT [%s]: t_enc=%.4f s  t_apply=%.4f s  total=%.4f s  Dmax=%d  k=%d\n",
                        COMP_LABEL, string(kind), stats_comp.time, stats_apply.time, stats_total.time, d_zt, ksig)
            else
                @printf("  %s+zT [%s]: t_enc=%.4f s  t_apply=%.4f s  total=%.4f s  Dmax=%d\n",
                        COMP_LABEL, string(kind), stats_comp.time, stats_apply.time, stats_total.time, d_zt)
            end

            save_state(series)

            if stats_total.time > TIME_TO_STOP
                skip[k1] = true
                skip[k2] = true
                println("  $(string(kind)) past time_to_stop; skipping larger n for this signal.")
            end
        end
    end

    save_state(series)
    println("\nArtifact: ", ART_PATH)
    return nothing
end

run()
