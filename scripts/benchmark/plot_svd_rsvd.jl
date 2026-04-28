# scripts/benchmark/plot_svd_rsvd.jl
#
# Fig A plotter. Single twin-axis SVG: mean runtime (left, solid) and peak
# memory (right, dashed) for ITensors.svd vs QILaplace.RSVD.rsvd, from
# scripts/benchmark/results/svd_rsvd_itensor.jld2.

using CairoMakie
using FileIO
using JLD2
using LaTeXStrings

include(joinpath(@__DIR__, "common.jl"))
include(joinpath(@__DIR__, "plot_utils.jl"))

const ART_PATH = joinpath(RESULTS_DIR, "svd_rsvd_itensor.jld2")
const OUT_PATH = joinpath(FIGURES_DIR, "svd_rsvd_time_mem.svg")

function run()
    data = load_results(ART_PATH)
    isempty(data) && error("No benchmark data at $ART_PATH. Run svd_rsvd_itensor.jl first.")

    svd_s  = from_dict(data["svd"])
    rsvd_s = from_dict(data["rsvd"])

    ns_svd  = sorted_ns(svd_s)
    ns_rsvd = sorted_ns(rsvd_s)

    time_svd  = series_vector(svd_s,  :time, ns_svd)
    time_rsvd = series_vector(rsvd_s, :time, ns_rsvd)
    mem_svd   = series_vector(svd_s,  :mem,  ns_svd)
    mem_rsvd  = series_vector(rsvd_s, :mem,  ns_rsvd)

    twin_axis_time_memory_plot(
        ns_svd, time_svd, ns_rsvd, time_rsvd,
        COLOR_SVD, COLOR_RSVD,
        L"n \, (\mathrm{qubits})",
        L"\mathrm{Time\ (s)}",
        L"\mathrm{Memory\ (MB)}",
        ns_svd, mem_svd, ns_rsvd, mem_rsvd;
        label_time_a = "SVD",
        label_time_b = "RSVD",
        outpath = OUT_PATH,
    )
    return nothing
end

run()
