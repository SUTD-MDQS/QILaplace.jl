# scripts/benchmark/plot_utils.jl
#
# Shared CairoMakie helpers used by the plot_*.jl scripts. Include this AFTER
# `using CairoMakie`/`using LaTeXStrings` and AFTER `include("common.jl")`.

"""
    twin_axis_plot(xs_a, ys_a, xs_b, ys_b, col_a, col_b,
                   xlabel, ylabel, label_a, label_b,
                   ylabel_secondary,
                   gcx_a, gcy_a, gcx_b, gcy_b,
                   title;
                   primary_log=true, secondary_log=true, outpath)

Draw two primary series (a, b) plus two dashed secondary series (gc_a, gc_b)
against the same x-axis using a right-side twin axis, then save to
`outpath`.
"""
function twin_axis_plot(
    xs_a, ys_a, xs_b, ys_b,
    col_a::AbstractString, col_b::AbstractString,
    xlabel::AbstractString, ylabel::AbstractString,
    label_a::AbstractString, label_b::AbstractString,
    ylabel_secondary::AbstractString,
    gcx_a, gcy_a, gcx_b, gcy_b,
    title::AbstractString;
    primary_log::Bool   = true,
    secondary_log::Bool = true,
    outpath::AbstractString,
)
    fig = Figure(; size = (1100, 700))

    ax_main = Axis(fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel,
        title  = title,
        yscale = primary_log ? log10 : identity,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        xlabelsize = 22, xticklabelsize = 22,
        ylabelsize = 22, yticklabelsize = 22,
        titlesize  = 20,
    )
    ax_main.yminorticks = IntervalsBetween(10)

    l_a = lines!(ax_main, xs_a, ys_a; color = col_a, linewidth = 2)
    s_a = scatter!(ax_main, xs_a, ys_a; color = col_a, marker = :circle, markersize = 12)
    l_b = lines!(ax_main, xs_b, ys_b; color = col_b, linewidth = 2)
    s_b = scatter!(ax_main, xs_b, ys_b; color = col_b, marker = :rect,   markersize = 12)

    ax_sec = Axis(fig[1, 1];
        xlabelvisible = false,
        xticksvisible = false,
        xticklabelsvisible = false,
        ylabel = ylabel_secondary,
        yscale = secondary_log ? log10 : identity,
        yaxisposition = :right,
        backgroundcolor = (:white, 0.0),
        xgridvisible = false, ygridvisible = false,
        xminorgridvisible = false, yminorgridvisible = false,
        ylabelsize = 22, yticklabelsize = 22,
    )

    allx = vcat(collect(xs_a), collect(xs_b), collect(gcx_a), collect(gcx_b))
    if !isempty(allx)
        xmin, xmax = extrema(allx)
        xlims!(ax_main, xmin, xmax)
        xlims!(ax_sec,  xmin, xmax)
    end

    l_gc_a = nothing
    if !isempty(gcx_a)
        valid = findall(y -> !isnan(y) && y > 0, gcy_a)
        if !isempty(valid)
            l_gc_a = lines!(ax_sec, gcx_a[valid], gcy_a[valid];
                             color = col_a, linewidth = 1.5, linestyle = :dash)
        end
    end
    l_gc_b = nothing
    if !isempty(gcx_b)
        valid = findall(y -> !isnan(y) && y > 0, gcy_b)
        if !isempty(valid)
            l_gc_b = lines!(ax_sec, gcx_b[valid], gcy_b[valid];
                             color = col_b, linewidth = 1.5, linestyle = :dash)
        end
    end

    artists = Vector{Vector{Any}}()
    labels  = String[]
    push!(artists, [l_a, s_a]); push!(labels, label_a)
    push!(artists, [l_b, s_b]); push!(labels, label_b)
    if l_gc_a !== nothing
        push!(artists, [l_gc_a]); push!(labels, "$(label_a) — $(ylabel_secondary)")
    end
    if l_gc_b !== nothing
        push!(artists, [l_gc_b]); push!(labels, "$(label_b) — $(ylabel_secondary)")
    end

    axislegend(ax_main, artists, labels;
               position = :lt, framevisible = true, labelsize = 18,
               patchsize = (18, 8), nbanks = 2)

    mkpath(dirname(outpath))
    save(outpath, fig)
    println("Saved: ", outpath)
    return fig
end

"""
    twin_axis_time_memory_plot(xs_a, ys_time_a, xs_b, ys_time_b, col_a, col_b,
                               xlabel, ylabel_time, ylabel_mem,
                               memx_a, mem_y_a, memx_b, mem_y_b; outpath, ...)

Single figure: **time** (left, solid lines + markers) and **memory** (right,
dashed lines) for two methods (e.g. SVD vs RSVD). No title. Both y-axes use
log scale by default.
"""
function twin_axis_time_memory_plot(
    xs_a, ys_time_a, xs_b, ys_time_b,
    col_a::AbstractString, col_b::AbstractString,
    xlabel, ylabel_time, ylabel_mem,
    memx_a, mem_y_a, memx_b, mem_y_b;
    label_time_a::AbstractString = "SVD",
    label_time_b::AbstractString = "RSVD",
    primary_log::Bool   = true,
    secondary_log::Bool = true,
    outpath::AbstractString,
)
    fig = Figure(; size = (1100, 700))

    ax_main = Axis(fig[1, 1];
        xlabel = xlabel,
        ylabel = ylabel_time,
        yscale = primary_log ? log10 : identity,
        xgridvisible = true, ygridvisible = true,
        xminorgridvisible = true, yminorgridvisible = true,
        xlabelsize = 22, xticklabelsize = 22,
        ylabelsize = 22, yticklabelsize = 22,
    )
    ax_main.yminorticks = IntervalsBetween(10)

    l_ta = lines!(ax_main, xs_a, ys_time_a; color = col_a, linewidth = 2)
    s_ta = scatter!(ax_main, xs_a, ys_time_a; color = col_a, marker = :circle, markersize = 12)
    l_tb = lines!(ax_main, xs_b, ys_time_b; color = col_b, linewidth = 2)
    s_tb = scatter!(ax_main, xs_b, ys_time_b; color = col_b, marker = :rect, markersize = 12)

    ax_sec = Axis(fig[1, 1];
        xlabelvisible = false,
        xticksvisible = false,
        xticklabelsvisible = false,
        ylabel = ylabel_mem,
        yscale = secondary_log ? log10 : identity,
        yaxisposition = :right,
        backgroundcolor = (:white, 0.0),
        xgridvisible = false, ygridvisible = false,
        xminorgridvisible = false, yminorgridvisible = false,
        ylabelsize = 22, yticklabelsize = 22,
    )
    ax_sec.yminorticks = IntervalsBetween(10)

    allx = vcat(collect(xs_a), collect(xs_b), collect(memx_a), collect(memx_b))
    if !isempty(allx)
        xmin, xmax = extrema(allx)
        xlims!(ax_main, xmin, xmax)
        xlims!(ax_sec, xmin, xmax)
    end

    l_ma = nothing
    if !isempty(memx_a)
        valid = findall(y -> !isnan(y) && y > 0, mem_y_a)
        if !isempty(valid)
            l_ma = lines!(ax_sec, memx_a[valid], mem_y_a[valid];
                          color = col_a, linewidth = 1.5, linestyle = :dash)
        end
    end
    l_mb = nothing
    if !isempty(memx_b)
        valid = findall(y -> !isnan(y) && y > 0, mem_y_b)
        if !isempty(valid)
            l_mb = lines!(ax_sec, memx_b[valid], mem_y_b[valid];
                          color = col_b, linewidth = 1.5, linestyle = :dash)
        end
    end

    # axislegend fills row-major with nbanks=2: entries 1–2 on row 1, 3–4 on row 2.
    # Order [SVD time, SVD mem, RSVD time, RSVD mem] ⇒ col1 = times, col2 = memory.
    artists = Vector{Vector{Any}}()
    labels  = String[]
    push!(artists, [l_ta, s_ta]); push!(labels, "$(label_time_a) time")
    if l_ma !== nothing
        push!(artists, [l_ma]); push!(labels, "$(label_time_a) memory")
    end
    push!(artists, [l_tb, s_tb]); push!(labels, "$(label_time_b) time")
    if l_mb !== nothing
        push!(artists, [l_mb]); push!(labels, "$(label_time_b) memory")
    end

    n_banks = length(artists) == 4 ? 2 : 1
    axislegend(ax_main, artists, labels;
               position = :lt, framevisible = true, labelsize = 18,
               patchsize = (18, 8), nbanks = n_banks)

    mkpath(dirname(outpath))
    save(outpath, fig)
    println("Saved: ", outpath)
    return fig
end
