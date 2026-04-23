# scripts/benchmark/plot_qft_vs_fftw.jl
#
# Fig D plotter. One log-y runtime panel with requested styling:
# - solid  : QFT only
# - dashdot : RSVD + QFT
# - black diamonds: FFTW(random) single baseline

using CairoMakie
using FileIO
using JLD2
using LaTeXStrings

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "qft_vs_fftw.jld2")
const OUT_PATH = joinpath(FIGURES_DIR, "qft_vs_fftw.svg")

const SIGNAL_KINDS = [:random, :sin, :sine20, :sin_cusp]

const LABEL_SHORT = Dict(
    :random   => "Random data",
    :sin      => "1 sine",
    :sine20   => "20 sines",
    :sin_cusp => "1 sine with cusp",
)

function run()
    data = load_results(ART_PATH)
    isempty(data) && error("No data at $ART_PATH. Run qft_vs_fftw.jl first.")

    fig = Figure(; size = (1000, 640))
    ax = Axis(fig[1, 1];
        xlabel = L"n \, (\mathrm{qubits})",
        ylabel = L"\mathrm{Time\ (s)}",
        yscale = log10,
        yminorticksvisible = true,
        xticks = 2:2:30,
        xlabelsize = 22, xticklabelsize = 20,
        ylabelsize = 22, yticklabelsize = 20,
        titlesize  = 20,
    )
    ylims!(ax, 1e-4, 1e1)
    ax.yminorticks = IntervalsBetween(10)

    artists = Vector{Vector{Any}}()
    labels  = Any[]

    for kind in SIGNAL_KINDS
        color  = get(SIGNAL_COLORS,  kind, "#000000")
        marker = get(SIGNAL_MARKERS, kind, :circle)

        key_qft_only = "qft_only_$(kind)"
        key_rsvd_qft = "rsvd_qft_$(kind)"
        haskey(data, key_qft_only) || continue
        haskey(data, key_rsvd_qft) || continue

        ser_qft_only = from_dict(data[key_qft_only])
        ser_rsvd_qft = from_dict(data[key_rsvd_qft])

        ns_qft_only = sorted_ns(ser_qft_only)
        ns_rsvd_qft = sorted_ns(ser_rsvd_qft)
        isempty(ns_qft_only) && isempty(ns_rsvd_qft) && continue

        t_qft_only = series_vector(ser_qft_only, :time, ns_qft_only)
        t_rsvd_qft = series_vector(ser_rsvd_qft, :time, ns_rsvd_qft)

        l_qft = lines!(ax, ns_qft_only, t_qft_only; color = color, linewidth = 2.5, linestyle = :solid)
        s_qft = scatter!(ax, ns_qft_only, t_qft_only; color = color, marker = marker, markersize = 12)
        l_rsvd = lines!(ax, ns_rsvd_qft, t_rsvd_qft; color = color, linewidth = 2.5, linestyle = :dashdot)

        push!(artists, [l_qft, s_qft]);   push!(labels, LaTeXString("\\mathrm{" * replace(LABEL_SHORT[kind], " " => "\\ ") * "\\,|\\,QFT\\ only}"))
        push!(artists, [l_rsvd]); push!(labels, LaTeXString("\\mathrm{" * replace(LABEL_SHORT[kind], " " => "\\ ") * "\\,|\\,RSVD{+}QFT}"))
    end

    if haskey(data, "fftw_random")
        ser_fft = from_dict(data["fftw_random"])
        ns_fft = sorted_ns(ser_fft)
        t_fft = series_vector(ser_fft, :time, ns_fft)
        l_fft = lines!(ax, ns_fft, t_fft; color = COLOR_FFTW, linewidth = 2.4, linestyle = :solid)
        s_fft = scatter!(ax, ns_fft, t_fft; color = COLOR_FFTW, marker = :diamond, markersize = 13)
        push!(artists, [l_fft, s_fft])
        push!(labels, L"\mathrm{FFTW\ (random\ only)}")
    end

    axislegend(ax, artists, labels;
        position = :lt,
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
