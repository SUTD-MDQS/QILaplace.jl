# QILaplace Benchmarking Suite

Standalone benchmark scripts that generate the data and figures embedded in
[`docs/src/benchmarking.md`](../../docs/src/benchmarking.md). These scripts are
intentionally **not** executed by `Literate.jl` during the docs build; instead
they are run locally on a reference machine (MacBook Pro M2 in our case) and
their outputs (JLD2 data + SVG figures) are committed to the repository.

## Folder Layout

```
scripts/benchmark/
  common.jl              # shared helpers (metadata, signal kinds, save/load, plot style)
  svd_rsvd_itensor.jl    # Fig A runner  : SVD vs RSVD on a random ITensor
  tt_decomp.jl           # Fig B runner  : signal_mps(:svd) vs signal_mps(:rsvd)
  mpo_bond_dim.jl        # Fig C runner  : MPO max bond dim for QFT / DT / zT
  qft_vs_fftw.jl         # Fig D runner  : end-to-end QFT vs FFTW bfft
  zt_full_runtime.jl     # Fig E runner  : signal_ztmps + apply(dt_mpo, .) per signal kind
  plot_bench_signals.jl  # Fig 1 plotter : time-domain view of every benchmark signal kind
  plot_svd_rsvd.jl       # Fig A plotter
  plot_tt_decomp.jl      # Fig B plotter
  plot_mpo_bond_dim.jl   # Fig C plotter
  plot_qft_vs_fftw.jl    # Fig D plotter
  plot_zt_runtime.jl     # Fig E plotter
  results/               # JLD2 artifacts from each runner
```

Plots are written to `docs/src/assets/benchmarking/`.

## Running

From the repository root:

```bash
julia --project=scripts/benchmark -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=scripts/benchmark scripts/benchmark/svd_rsvd_itensor.jl
julia --project=scripts/benchmark scripts/benchmark/tt_decomp.jl
julia --project=scripts/benchmark scripts/benchmark/mpo_bond_dim.jl
julia --project=scripts/benchmark scripts/benchmark/qft_vs_fftw.jl
julia --project=scripts/benchmark scripts/benchmark/zt_full_runtime.jl
```

Each runner writes its JLD2 file into `scripts/benchmark/results/` and can be
re-invoked to append/extend data — a `time_to_stop` guard skips further
data points once a method crosses a hard time budget, so the same script can be
run on larger sweeps without redoing the work.

Then generate the figures:

```bash
julia --project=scripts/benchmark scripts/benchmark/plot_bench_signals.jl
julia --project=scripts/benchmark scripts/benchmark/plot_svd_rsvd.jl
julia --project=scripts/benchmark scripts/benchmark/plot_tt_decomp.jl
julia --project=scripts/benchmark scripts/benchmark/plot_mpo_bond_dim.jl
julia --project=scripts/benchmark scripts/benchmark/plot_qft_vs_fftw.jl
julia --project=scripts/benchmark scripts/benchmark/plot_zt_runtime.jl
```

## Data Schema (per JLD2 file)

Every runner stores a uniform schema:

```julia
meta = Dict(
    "julia_version" => string(VERSION),
    "blas_threads"  => BLAS.get_num_threads(),
    "cpu_brand"     => Sys.cpu_info()[1].model,
    "date"          => string(now()),
    "params"        => Dict("cutoff" => ..., "k" => ..., ...),
)
# plus one or more result series keyed by a descriptive symbol/string
```

For runtime scripts the series carry `("n", "time", "gctime", "mem", "allocs", "maxbond")`.
For the zT full runtime script each signal kind stores `("n", "t_input_mps", "t_apply")`.
