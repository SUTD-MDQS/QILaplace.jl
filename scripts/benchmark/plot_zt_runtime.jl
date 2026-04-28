# scripts/benchmark/plot_zt_runtime.jl
#
# Fig E plotter. Log-y runtime: solid = apply(W, ψ) only, dashed = encode + apply.
# x-axis: output qubit count m = 2n. Legend layout matches plot_qft_vs_fftw.jl
# (nbanks = 2, concise math-style labels).

using CairoMakie
using FileIO
using JLD2
using LaTeXStrings

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "zt_full_runtime.jld2")
const OUT_PATH = joinpath(FIGURES_DIR, "zt_full_runtime.svg")

const SIGNAL_KINDS = [:sin, :multi_sin, :multi_sin_exp, :abs_cos_power_p8, :random]

# Short names (spaces → \  in \mathrm{...}) — same spirit as LABEL_SHORT in plot_qft_vs_fftw.jl
const LABEL_SHORT = Dict(
    :sin              => "1 sine",
    :multi_sin_exp    => "Damped sine",
    :abs_cos_power_p8 => "Abs cos 0.8",
    :random           => "Random",
)

function zt_legend_label(short::AbstractString, suffix::AbstractString)
    esc = replace(short, " " => "\\ ")
    return LaTeXString("\\mathrm{" * esc * "\\,|\\," * suffix * "}")
end

# Decode the ZTSeries schema written by zt_full_runtime.jl
struct ZTPoints
    ns::Vector{Int}
    t_input_mps::Vector{Float64}
    t_apply::Vector{Float64}
end

function load_zt(data::AbstractDict, kind::Symbol)
    label = String(kind)
    haskey(data, label) || return ZTPoints(Int[], Float64[], Float64[])
    d = data[label]
    ns = sort(collect(keys(d["t_apply"])))
    return ZTPoints(
        ns,
        Float64[d["t_input_mps"][n] for n in ns],
        Float64[d["t_apply"][n]     for n in ns],
    )
end

function run()
    data = load_results(ART_PATH)
    isempty(data) && error("No data at $ART_PATH. Run zt_full_runtime.jl first.")

    xtick_vals = 2*collect(2:4:30)
    xtick_labs = string.(xtick_vals)

    ylo, yhi = 1e-4, 1e2 * 2.0
    k_lo, k_hi = -4, 2
    ytick_vals = [10.0^k for k in k_lo:2:k_hi]
    ytick_labs = [LaTeXString("10^{$k}") for k in k_lo:2:k_hi]

    fig = Figure(; size = (1000, 640))
    ax = Axis(fig[1, 1];
        xlabel = L"m \, (\mathrm{output\ qubits})",
        xticks = (xtick_vals, xtick_labs),
        ylabel = L"\mathrm{Time\ (s)}",
        yscale = log10,
        yticks = (ytick_vals, ytick_labs),
        yminorticksvisible = true,
        xlabelsize = 22, xticklabelsize = 20,
        ylabelsize = 22, yticklabelsize = 20,
        titlesize  = 20,
    )
    ylims!(ax, ylo, yhi)
    ax.yminorticks = IntervalsBetween(10)

    artists = Vector{Vector{Any}}()
    labels  = Any[]

    for kind in SIGNAL_KINDS
        pts = load_zt(data, kind)
        isempty(pts.ns) && continue

        m = 2 .* pts.ns
        t_core = pts.t_apply
        t_full = pts.t_input_mps .+ pts.t_apply
        color  = get(SIGNAL_COLORS,  kind, "#000000")
        marker = get(SIGNAL_MARKERS, kind, :circle)

        l_core = lines!(ax, m, t_core; color = color, linewidth = 2.5, linestyle = :solid)
        s_core = scatter!(ax, m, t_core; color = color, marker = marker, markersize = 12)
        l_full = lines!(ax, m, t_full; color = color, linewidth = 2.5, linestyle = :dash, alpha = 0.75)

        short = LABEL_SHORT[kind]
        # Row-major nbanks=2: [apply, full] per kind ⇒ col1 = apply, col2 = enc+zT
        push!(artists, [l_core, s_core])
        push!(labels, zt_legend_label(short, "apply"))
        push!(artists, [l_full])
        push!(labels, zt_legend_label(short, "enc{+}zT"))
    end

    axislegend(ax, artists, labels;
        position = :rb,
        framevisible = true,
        labelsize = 13,
        patchsize = (24, 10),
        nbanks = 2,
        margin = (10, 10, 10, 10),
    )

    save(OUT_PATH, fig)
    println("Saved: ", OUT_PATH)
    return nothing
end

run()
