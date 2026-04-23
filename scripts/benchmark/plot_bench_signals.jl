# scripts/benchmark/plot_bench_signals.jl
#
# Figure 1 plotter. Two-panel figure showing the benchmark signal families
# in the time domain at n = 10 (N = 1024 samples):
# - Top:    signals used by qft_vs_fftw.jl   (SIGNAL_KINDS there)
# - Bottom: signals used by zt_full_runtime.jl (SIGNAL_KINDS there)
#
# Output: docs/src/assets/benchmarking/bench_signals.svg

using CairoMakie
using LaTeXStrings

include(joinpath(@__DIR__, "common.jl"))

const OUT_PATH = joinpath(FIGURES_DIR, "bench_signals.svg")

const N_QUBITS = 10 # N = 1024

const QFT_KINDS = [:sin, :sine20, :sin_cusp, :random]
const ZT_KINDS  = [:sin, :multi_sin_exp, :abs_cos_power_p8, :random]

const LABELS = Dict(
    :sin              => "1 sine",
    :sine20           => "20 sines",
    :sin_cusp         => "1 sine + cusps",
    :multi_sin        => "Multi-sine",
    :multi_sin_exp    => "Damped multi-sine",
    :abs_cos_power_p8 => L"|\cos|^{0.8}",
    :random           => "Random",
)

function plot_panel!(ax, kinds::Vector{Symbol}, n::Int; legend_position::Symbol = :rt)
    artists = Any[]
    labels  = Any[]
    N = 2^n
    xs = 0:(N - 1)
    # Draw random first so structured signals stay visible on top.
    draw_kinds = Symbol[]
    if :random in kinds
        push!(draw_kinds, :random)
    end
    append!(draw_kinds, filter(!=(:random), kinds))

    for kind in draw_kinds
        y = make_signal(kind, n)
        color = get(SIGNAL_COLORS, kind, "#000000")
        if kind == :random
            l = lines!(ax, xs, y; color = (color, 0.25), linewidth = 1.0)
        else
            l = lines!(ax, xs, y; color = color, linewidth = 1.6)
        end
        push!(artists, l)
        push!(labels, LABELS[kind])
    end
    axislegend(ax, artists, labels;
        position   = legend_position,
        framevisible = true,
        labelsize  = 14,
        patchsize  = (20, 8),
        nbanks     = 2,
        margin     = (6, 6, 6, 6),
    )
    return nothing
end

function run()
    N = 2^N_QUBITS
    fig = Figure(; size = (1120, 560), figure_padding = (8, 12, 8, 12))

    ax_top = Axis(fig[1, 1];
        title  = "QFT benchmark signals",
        xlabel = L"j \, (\mathrm{sample\ index})",
        ylabel = L"x_j",
        xticks = 0:(N >> 3):N,
        xlabelsize = 18, xticklabelsize = 14,
        ylabelsize = 18, yticklabelsize = 14,
        titlesize  = 18,
    )

    ax_bot = Axis(fig[2, 1];
        title  = "Discrete Laplace (zT) benchmark signals",
        xlabel = L"j \, (\mathrm{sample\ index})",
        ylabel = L"x_j",
        xticks = 0:(N >> 3):N,
        xlabelsize = 18, xticklabelsize = 14,
        ylabelsize = 18, yticklabelsize = 14,
        titlesize  = 18,
    )

    plot_panel!(ax_top, QFT_KINDS, N_QUBITS; legend_position = :lb)
    plot_panel!(ax_bot, ZT_KINDS,  N_QUBITS; legend_position = :rb)

    ylims!(ax_top, -2, 2)
    ylims!(ax_bot, -2, 2)
    xlims!(ax_top, 0, N - 1)
    xlims!(ax_bot, 0, N - 1)

    rowgap!(fig.layout, 12)

    save(OUT_PATH, fig)
    println("Saved: ", OUT_PATH)
    return nothing
end

run()
