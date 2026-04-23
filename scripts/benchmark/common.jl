# scripts/benchmark/common.jl
# Shared helpers for the QILaplace benchmarking suite: machine metadata
# capture, signal generators (including kinds not exposed by the public
# `generate_signal` API), a uniform save/load schema for the per-figure
# JLD2 artifacts, and a consistent CairoMakie style palette.

using BenchmarkTools
using Dates
using FileIO
using JLD2
using LinearAlgebra
using Printf
using Random
using Statistics

using QILaplace

const RESULTS_DIR = joinpath(@__DIR__, "results")
const FIGURES_DIR = normpath(joinpath(@__DIR__, "..", "..", "docs", "src", "assets", "benchmarking"))

mkpath(RESULTS_DIR)
mkpath(FIGURES_DIR)

# ------------------------------------------------------------------
# Machine / software metadata
# ------------------------------------------------------------------

function machine_meta(; params::Dict = Dict{String, Any}())
    cpu = try
        Sys.cpu_info()[1].model
    catch
        "unknown"
    end
    return Dict{String, Any}(
        "julia_version" => string(VERSION),
        "blas_threads"  => BLAS.get_num_threads(),
        "num_threads"   => Threads.nthreads(),
        "cpu_brand"     => cpu,
        "sys_total_mem_gb" => round(Sys.total_memory() / 2^30; digits = 2),
        "date"          => string(now()),
        "params"        => Dict{String, Any}(String(k) => v for (k, v) in params),
    )
end

# ------------------------------------------------------------------
# Signal generators
# ------------------------------------------------------------------
# Wraps QILaplace.generate_signal for benchmark signal families.

"""
    make_signal(kind, n; seed=1234)

Return a length-`2^n` real vector of the requested shape. Supports:
- `:sin`       — single-frequency sine via `generate_signal`
- `:multi_sin`, `:multi_sin_exp`, `:abs_cos_power_p8` — generated via `generate_signal` (deterministic parameterisation for multi-signal variants)
- `:sine20`    — multi-sine via `generate_signal` with 20 radial frequencies `2πk/N`
- `:sin_cusp`  — one cosine on `[0,1]` plus four terms `a e^{-b|x-c|}` (sharp cusps)
- `:random`    — Gaussian noise, reproducible via `seed`
"""
# Cusp benchmark: cos(2πx) + Σ a_i exp(-b_i |x - c_i|), x ∈ [0,1]
const _CUSP_A = Float64[0.6, -0.4, 0.5, 0.35]
const _CUSP_B = Float64[80.0, 120.0, 90.0, 100.0]
const _CUSP_C = Float64[0.2, 0.45, 0.62, 0.85]

function make_signal(kind::Symbol, n::Int; seed::Int = 1234)
    N = 2^n
    if kind === :sin
        return generate_signal(n; kind = :sin, dt = 1.0, freq = 2π * 2 / N)
    elseif kind === :multi_sin
        return generate_signal(n; kind = :multi_sin, dt = 5.0 / N)
    elseif kind === :multi_sin_exp
        return generate_signal(n; kind = :multi_sin_exp, dt = 5.0 / N, ω_scale = 150.0)
    elseif kind === :abs_cos_power_p8
        return generate_signal(n; kind = :abs_cos_power_p8, dt = 5.0 / N)
    elseif kind === :sine20
        freqs = [2π * k / N for k in 1:20]
        return generate_signal(n; kind = :sin, dt = 1.0, freq = freqs)
    elseif kind === :sin_cusp
        return [
            let x = N == 1 ? 0.0 : j / (N - 1)
                cos(2π * x) +
                sum(ai * exp(-bi * abs(x - ci)) for (ai, bi, ci) in zip(_CUSP_A, _CUSP_B, _CUSP_C))
            end for j in 0:(N - 1)
        ]
    elseif kind === :random
        rng = Xoshiro(seed + n)
        return randn(rng, N)
    else
        throw(ArgumentError("make_signal: unknown kind $kind"))
    end
end

# ------------------------------------------------------------------
# Benchmark stats extraction (BenchmarkTools.Trial)
# ------------------------------------------------------------------

struct RunStats
    time::Float64     # seconds
    gctime::Float64   # seconds
    mem::Float64      # MB
    allocs::Int       # count
end

function stats_from_trial(t::BenchmarkTools.Trial)
    m = mean(t)
    return RunStats(m.time / 1e9, m.gctime / 1e9, m.memory / 1e6, Int(m.allocs))
end

"""
    elapsed_stats(thunk) -> RunStats

Fallback for expensive runs where a full `@benchmark` sample set is too costly.
Uses `@timed` for a single run (after a prior warmup run by the caller).
"""
function elapsed_stats(thunk)
    result = @timed thunk()
    return RunStats(result.time, result.gctime, result.bytes / 1e6, Int(result.compile_time > 0 ? 0 : 0))
end

# ------------------------------------------------------------------
# Series containers (what each runner writes to JLD2)
# ------------------------------------------------------------------

"""
    Series()

Empty container for `(n, time, gctime, mem, allocs, maxbond)` benchmark series.
Dictionaries are keyed by `n` so the scripts can resume incrementally.
"""
mutable struct Series
    time::Dict{Int, Float64}
    gctime::Dict{Int, Float64}
    mem::Dict{Int, Float64}
    allocs::Dict{Int, Int}
    maxbond::Dict{Int, Int}
end

Series() = Series(
    Dict{Int, Float64}(),
    Dict{Int, Float64}(),
    Dict{Int, Float64}(),
    Dict{Int, Int}(),
    Dict{Int, Int}(),
)

function update!(s::Series, n::Int, stats::RunStats; maxbond::Int = 0)
    s.time[n] = stats.time
    s.gctime[n] = stats.gctime
    s.mem[n] = stats.mem
    s.allocs[n] = stats.allocs
    s.maxbond[n] = maxbond
    return s
end

has_point(s::Series, n::Int) = haskey(s.time, n)

to_dict(s::Series) = Dict{String, Any}(
    "time"    => s.time,
    "gctime"  => s.gctime,
    "mem"     => s.mem,
    "allocs"  => s.allocs,
    "maxbond" => s.maxbond,
)

function from_dict(d::AbstractDict)
    s = Series()
    haskey(d, "time")    && merge!(s.time,    d["time"])
    haskey(d, "gctime")  && merge!(s.gctime,  d["gctime"])
    haskey(d, "mem")     && merge!(s.mem,     d["mem"])
    haskey(d, "allocs")  && merge!(s.allocs,  d["allocs"])
    haskey(d, "maxbond") && merge!(s.maxbond, d["maxbond"])
    return s
end

sorted_ns(s::Series) = sort(collect(keys(s.time)))

function series_vector(s::Series, field::Symbol, ns::AbstractVector{<:Integer})
    d = getfield(s, field)
    T = valtype(d)
    return T[get(d, n, T === Int ? 0 : NaN) for n in ns]
end

# ------------------------------------------------------------------
# JLD2 I/O
# ------------------------------------------------------------------

"""
    save_results(path; meta, series_dict)

Save a uniform benchmark artifact. `series_dict` maps a string label (e.g.
`"svd"`, `"rsvd"`, `"qft"`) to either a `Series` or an arbitrary Dict.
"""
function save_results(path::AbstractString; meta::AbstractDict, series::AbstractDict)
    mkpath(dirname(path))
    jldopen(path, "w") do f
        f["meta"] = Dict{String, Any}(String(k) => v for (k, v) in meta)
        for (label, ser) in series
            d = ser isa Series ? to_dict(ser) : ser
            f[String(label)] = d
        end
    end
    return path
end

function load_results(path::AbstractString)
    return isfile(path) ? load(path) : Dict{String, Any}()
end

function load_series(path::AbstractString, label::AbstractString)
    data = load_results(path)
    return haskey(data, label) ? from_dict(data[label]) : Series()
end

# ------------------------------------------------------------------
# Max bond-dim helpers
# ------------------------------------------------------------------

maxbond_mps(ψ::SignalMPS) = isempty(ψ.bonds) ? 1 : maximum(dim, ψ.bonds)

function maxbond_mps(ψ::zTMPS)
    ds = Int[]
    append!(ds, dim.(ψ.bonds_main))
    append!(ds, dim.(ψ.bonds_copy))
    return isempty(ds) ? 1 : maximum(ds)
end

function maxbond_mpo(W)
    # Works for both SingleSiteMPO and PairedSiteMPO by inspecting fields.
    ds = Int[]
    if hasproperty(W, :bonds)
        append!(ds, dim.(W.bonds))
    end
    if hasproperty(W, :bonds_main)
        append!(ds, dim.(W.bonds_main))
    end
    if hasproperty(W, :bonds_copy)
        append!(ds, dim.(W.bonds_copy))
    end
    return isempty(ds) ? 1 : maximum(ds)
end

# ------------------------------------------------------------------
# Logging
# ------------------------------------------------------------------

function banner(title::AbstractString)
    bar = "=" ^ max(length(title) + 4, 40)
    println()
    println(bar)
    println("  ", title)
    println(bar)
    return nothing
end

# ------------------------------------------------------------------
# Shared plotting palette (used by every plot_*.jl script)
# ------------------------------------------------------------------

const COLOR_SVD  = "#FF8C42"    # warm orange — slow / reference
const COLOR_RSVD = "#54B07E"    # green       — fast / randomised
const COLOR_FFTW = "#000000"    # black       — classical FFT baseline

const COLOR_QFT = "#D45E00"
const COLOR_DT  = "#0071B1"
const COLOR_ZT  = "#009E73"

const SIGNAL_COLORS = Dict(
    :sin              => "#0071B1",
    :multi_sin        => "#CC79A7",
    :multi_sin_exp    => "#D45E00",
    :abs_cos_power_p8 => "#009E73",
    :sine20           => "#D45E00",
    :sin_cusp         => "#009E73",
    :random           => "#9D4DC4",
)
const SIGNAL_MARKERS = Dict(
    :sin              => :circle,
    :multi_sin        => :diamond,
    :multi_sin_exp    => :rect,
    :abs_cos_power_p8 => :utriangle,
    :sine20           => :utriangle,
    :sin_cusp         => :rect,  # :rectangle is not reliably supported in CairoMakie/SVG
    :random           => :star5,
)
