# QILaplace.jl

<img src="logo.svg" alt="QILaplace logo" align="left" width="230"/>


QILaplace.jl is a Julia library for implementing **quantum-inspired signal transforms**—most notably the **Quantum Fourier Transform (QFT)** and **Discrete Laplace (aka z-Transform)** variants—using **tensor-network representations**. The package is built on **Matrix Product States** (MPS) and **Matrix Product Operators** (MPOs) and runs entirely on classical hardware via ITensors.jl.

The central idea is to leverage the representability of certain quantum circuits as MPOs and realize their action through efficient tensor-network contractions. This allows us to study, prototype, and apply QFT-like transforms at problem sizes far beyond what dense representations would permit—**without requiring a quantum computer**.

## Why tensor-network transforms?

### The classical state of the art

For discrete signal transforms, classical algorithms are extremely effective—but they implicitly assume **dense access** to all signal samples:

- **Discrete Fourier Transform (FFT)**\
  Runtime: **$O(N \ \log(N))$**\
  Requires storing and processing all $N$ samples explicitly.

- **Discrete Laplace / z-transform**\
  Naive implementations scale as **$O(N^2)$** and can suffer from [numerical instability](https://www.researchgate.net/publication/373077088_Fast_discrete_Laplace_transforms) at large $N$. However, several *fast classical alternatives* exist for specific evaluation settings:

  - **Fast evaluation on the positive real axis:** using specialized schemes, evaluation along the real axis can be performed in **effectively linear or near-linear time** with cost **$O(N + M)$**.
  - **Chirp-z transform (CZT):** reduces evaluation along specific contours in the complex plane to FFT-like convolutions, with cost **$O((N + M)\ \log(N + M))$**.
  - **Dense 2D evaluation grids:** typical scaling **$O(N·M)$**.

  These algorithms substantially outperform brute-force Laplace evaluation, but still fundamentally rely on *dense access to the signal* and do not exploit internal low-rank or entanglement structure.

However, many real-world signals are **structured and compressible**, with most of their "useful" information concentrated in a small number of dominant modes. Such signals admit low-rank representations and can be faithfully approximated with far fewer degrees of freedom than their ambient dimension suggests. QILaplace.jl adopts a tensor-network perspective:

- Signals of length $N$ are represented as **MPS with $\log_2(N)$ sites**
- Transform operators are constructed by efficiently contracting **logarithmic-depth MPO circuits**
- Approximation error, compression, and circuit structure are **explicitly controlled**

For signals admitting a bounded MPS bond dimension, this leads to a radically different scaling regime.


### What this package provides

QILaplace.jl implements three quantum-inspired discrete transforms:

- **Quantum Fourier Transform (QFT)**\
  A circuit-structured alternative to the FFT, realized as a finite-bond-dimension MPO.

- **Discrete Laplace Transform with real exponents (Damping Transform)**\
  Suitable for exponentially damped signals and decay analysis.

- **Full complex Discrete Laplace / z-Transform (zT)**\
  A tensor-network formulation of the classical z-transform, enabling pole–zero analysis directly from tensor states.


## Performance and scaling

The transform operators in QILaplace.jl correspond to **logarithmic-depth quantum circuits** represented as MPOs. Prior work of [Chen, Stoudenmire, and White](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.040318) has shown that QFT circuits admit an MPO representation with a **finite bond dimension** that does **not scale with system size** for fixed accuracy. We have observed the same behavior in the case of the discrete Laplace transform as well. 

Input signals are represented as MPS with bond dimension $\chi_s$, determined by the structure and complexity of the data. Structured or weakly correlated signals admit compact representations, while more complex signals may require larger bond dimensions.

Transforms are applied via MPO–MPS contraction. The computational cost scales as:

- Linear in the number of tensor sites: **$n = \log_2(N)$**
- Polynomial in the signal bond dimension $\chi_s$
- Polynomial in the circuit bond dimension $\chi_c$ (constant for fixed accuracy)

For fixed circuit accuracy, the effective runtime is **$O(\chi_s^2 \log(N))$**.

### Runtime comparison

| Transform                  | Best-known algorithm                        | Best-known classical runtime   | QILaplace runtime\* |
| -------------------------- | ------------------------------------------- | ------------------------------ | ------------------- |
| Discrete Fourier Transform | FFT (Cooley–Tukey)                          | $O(N \log(N))$                 | $O(\chi_s^2 \log(N))$        |
| Discrete Laplace (real)    | Chirp-z transform (CZT) / fast real-axis zT | $O(N + M)$ (best case)         | $O(\chi_s^2 \log(N))$        |
| z-Transform (complex)      | Chirp-z transform (CZT) on 2D grids         | $O(NM)$                        | $O(\chi_s^2 \log(N))$        |
| QFT (quantum circuit)      | Quantum QFT circuit                         | $O(\log^2(N))$†                | $O(\chi_s^2 \log(N))$        |

\* For signals admitting bounded MPS bond dimension.

† Idealized quantum circuit, ignoring state preparation, readout, and noise overheads.

In other words, **when the signal is compressible**, transform costs scale logarithmically with signal length. In this regime, tensor-network implementations can be competitive with—and in practice outperform—both dense classical algorithms and idealized quantum-circuit implementations for large data sizes, while running entirely on classical hardware. As usual, the bond dimension keeps score.


## Core features

### Signal generation and representation

- Signals are represented as MPS objects (e.g. `SignalMPS`, `zTMPS`)
- Dense arrays can be compressed using SVD or randomized SVD (rSVD)
- Built-in utilities generate oscillatory, damped, and structured test signals

### Transform construction

- QFT, damping-transform, and z-transform operators are constructed as MPOs
- Circuit-style MPO building blocks:
  - `SingleSiteMPO`
  - `PairedSiteMPO`
- High-level builders:
  - `build_qft_mpo`
  - `build_dt_mpo`
  - `build_zt_mpo`

### Transform application

- Generic interface: `transform(::AbstractMPS, ::AbstractMPO)`
- Outputs are produced in **bit-reversed order**, matching the natural QFT circuit layout

### Post-processing and analysis

- Pole and zero extraction directly from transformed tensor states
- Parameter estimation and spectral analysis tools (ongoing and planned)


## Authorship and maintenance

This repository supports the work reported in an upcoming arXiv manuscript.

- **Original author:** **Noufal Jaseem**\
  Developed the original implementation and core workflow used to generate the results reported in the manuscript.

- **Lead maintainer / major contributor:** **Gauthameshwar S.**\
  Curated and structured the repository, improved code quality and usability, and implemented targeted efficiency optimizations. The optimized implementation was used to update the runtime figures in the main manuscript.

- **Project supervision:** **Dario Poletti**\
  Provided scientific supervision and guidance for the project.


## Contributing

Yes—please don’t hesitate to contribute and help further develop this project.

This repository is maintained by a sleep-deprived PhD student who learned open-source development by starring popular Julia packages, reading forum threads, and watching YouTube tutorials. If parts of the codebase seem unpolished or unclear, please feel free to help me out.

Pull requests for refactors, documentation, performance improvements, or new ideas are very welcome. If you make the project better, you are a contributor 🫰🏻.


## Citation

If you use QILaplace.jl in published work, please cite the repository and the associated arXiv paper (link forthcoming).

Made with ❤️ using Julia, ITensors, and a healthy respect for bond-dimension growth.

---
