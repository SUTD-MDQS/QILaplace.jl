# scripts/benchmark/plot_mpo_bond_dim.jl
#
# Fig C plotter. Single-panel log-y plot of max MPO bond dimension versus
# qubit count, matching the left panel of scripts/Deprecated/plot_benchmark_time.jl.

using CairoMakie
using FileIO
using JLD2
using LaTeXStrings

include(joinpath(@__DIR__, "common.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "mpo_bond_dim.jld2")
const OUT_PATH = joinpath(FIGURES_DIR, "mpo_bond_dim.svg")

"""Minor ticks at 2·10^k … 9·10^k within [ylo, yhi] (standard log-scale grid)."""
function log_decade_minor_ticks(ylo, yhi, k_lo, k_hi)
    vals = Float64[]
    for k in k_lo:k_hi
        base = 10.0^k
        for m in 2:9
            v = m * base
            if ylo <= v <= yhi
                push!(vals, v)
            end
        end
    end
    return vals
end

function run()
    data = load_results(ART_PATH)
    isempty(data) && error("No data at $ART_PATH. Run mpo_bond_dim.jl first.")

    qft_s = from_dict(data["qft"])
    dt_s  = from_dict(data["dt"])
    zt_s  = from_dict(data["zt"])

    ns_qft = sorted_ns(qft_s)
    ns_dt  = sorted_ns(dt_s)
    ns_zt  = sorted_ns(zt_s)

    d_qft = [qft_s.maxbond[n] for n in ns_qft]
    d_dt  = [dt_s.maxbond[n]  for n in ns_dt]
    d_zt  = [zt_s.maxbond[n]  for n in ns_zt]

    # x-axis: input size n (qubits) for QFT, DT, and zT.
    x_qft = ns_qft
    x_dt  = ns_dt
    x_zt  = ns_zt

    all_d = vcat(d_qft, d_dt, d_zt)
    @assert all(all_d .> 0) "Bond dimensions must be positive for log scale."
    ymin = minimum(all_d)
    ylo  = 10.0 ^ floor(log10(ymin * 0.9))
    # Headroom above 10^2 so curves do not sit on the top spine.
    yhi  = 10.0^2 * 10^0.5

    # Major y ticks only at exact decades 10^k (no interpolated log ticks).
    k_hi = Int(ceil(log10(yhi)))
    while k_hi > Int(floor(log10(ylo))) && 10.0^k_hi > yhi
        k_hi -= 1
    end
    k_lo = Int(floor(log10(ylo)))
    ytick_vals = [10.0^k for k in k_lo:k_hi]
    ytick_labs = [latexstring("10^{$k}") for k in k_lo:k_hi]

    fig = Figure(; size = (700, 500))
    ax = Axis(fig[1, 1];
        xlabel = L"n \, (\mathrm{qubits})",
        ylabel = L"D_{\max}",
        yscale = log10,
        yticks = (ytick_vals, ytick_labs),
        yminorticks = log_decade_minor_ticks(ylo, yhi, k_lo, k_hi),
        xgridvisible = false, ygridvisible = false,
        yminorticksvisible = true,
        yticksize = 6, yminorticksize = 3,
        xlabelsize = 22, xticklabelsize = 20,
        ylabelsize = 22, yticklabelsize = 20,
        titlesize  = 20,
    )

    line_qft = lines!(ax, x_qft, d_qft; color = COLOR_QFT, linewidth = 2)
    scat_qft = scatter!(ax, x_qft, d_qft; color = COLOR_QFT, marker = :rect,      markersize = 13)
    line_dt  = lines!(ax, x_dt,  d_dt;  color = COLOR_DT,  linewidth = 2)
    scat_dt  = scatter!(ax, x_dt,  d_dt;  color = COLOR_DT,  marker = :circle,    markersize = 13)
    line_zt  = lines!(ax, x_zt,  d_zt;  color = COLOR_ZT,  linewidth = 2)
    scat_zt  = scatter!(ax, x_zt,  d_zt;  color = COLOR_ZT,  marker = :utriangle, markersize = 13)

    ylims!(ax, ylo, yhi)

    artists = [
        [line_qft, scat_qft],
        [line_dt,  scat_dt],
        [line_zt,  scat_zt],
    ]
    labels = [
        L"\mathrm{QFT}",
        L"\mathrm{DT}",
        L"\mathrm{zT}",
    ]
    axislegend(ax, artists, labels;
               position = :rb, framevisible = true,
               labelsize = 18, patchsize = (24, 10), nbanks = 1)

    save(OUT_PATH, fig)
    println("Saved: ", OUT_PATH)
    return nothing
end

run()
