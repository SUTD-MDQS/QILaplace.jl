<table align="center">
  <tr>
    <td style="padding-right: 16px;">
      <img src="logo.svg" alt="QILaplace logo" width="180">
    </td>
    <td>
      <h1>QILaplace.jl</h1>
    </td>
  </tr>
</table>

## Quantum-Inspired Laplace & Fourier transforms with tensor networks

QILaplace.jl is a Julia library for implementing quantum-inspired signal transforms—most notably the Quantum Fourier Transform (QFT) and Discrete Laplace (aka z-Transform) variants—using tensor-network representations. The package is built on Matrix Product States (MPS) and Matrix Product Operators (MPOs) and runs entirely on classical hardware via ITensors.jl.

The central idea is to leverage the structure of logarithmic-depth quantum circuits and realize them as efficient tensor networks. This allows us to study, prototype, and apply QFT-like transforms at problem sizes far beyond what dense representations would permit—without requiring a quantum computer.

This repository sits at the intersection of signal processing, quantum information, and numerical linear algebra. It is research-driven, performance-aware, and occasionally opinionated about bond dimensions.


## Why this package exists

Classical FFT-based pipelines are extremely effective, but they assume dense data representations. For large, structured, or compressible signals, tensor-network methods offer a complementary approach: explicit control over approximation error, compression, and circuit structure.

QILaplace.jl provides:
- Circuit-style constructions of QFT-, DT-, and z-transform operators as MPOs
- Efficient application of these operators to signal MPS representations
- Utilities for generating, compressing, and transforming signals
- Tools for pole and zero extraction directly from transformed tensor states

If you are curious about what happens when signal processing meets quantum circuit structure—this is the sandbox.


## Who is this for?

- Researchers in signal processing, computational physics, and quantum information
- Anyone experimenting with tensorized representations of large signals
- Julia users looking for an ITensors.jl-based QFT-style toolkit


## Core ideas and design

### Signal representation
- Signals are represented as MPS objects (e.g., `SignalMPS`, `zTMPS`)
- Dense arrays can be compressed using SVD or randomized SVD (rSVD)
- A convenience `generate_signal` function produces oscillatory and damped test signals

### Transform operators
- QFT, DT, and z-transform operators are constructed as MPOs
- Implemented using:
  - `SingleSiteMPO`
  - `PairedSiteMPO`
- Builder functions include:
  - `build_qft_mpo`
  - `build_dt_mpo`
  - `build_zt_mpo`

### Algorithms
- Generic application interface: `transform(::AbstractMPS, ::AbstractMPO)`
- Outputs are produced in **bit-reversed order**, matching the natural QFT circuit layout
- Planned and ongoing work includes:
    - bit-reversal back to physical ordering
    - pole and zero extraction
    - parameter estimation from transformed tensor states


## Performance and scaling

The transform operators in QILaplace.jl are implemented as Matrix Product Operators corresponding to logarithmic-depth quantum circuits. [Chen, Stoudenmire, and White (2023)](https://journals.aps.org/prxquantum/abstract/10.1103/PRXQuantum.4.040318) have shown that a QFT circuit MPO has a **finite bond dimension** $\chi_c$ that does **not scale with system size** for a fixed approximation accuracy. We have also shown similar bond dimension behavior for DT- and z-transform circuits (look for our paper on arxiv!).

Input signals are represented as Matrix Product States with bond dimension $\chi_s$, which depends on the structure and complexity of the data. Structured or weakly correlated signals admit compact representations, while more complex signals may require larger bond dimensions.

Transforms are applied via MPO–MPS contraction. The computational cost scales as:
- Linear in the number of tensor sites $n = \log_2 N$
- Polynomial in the signal bond dimension $\chi_s$
- Polynomial in the circuit bond dimension $\chi_c$ (which is constant)

For fixed circuit accuracy, the effective runtime is:

> **$O(\chi_s^2 \log(N))$**

In other words, when the input signal is well-represented by a finite-bond-dimension MPS, QFT-, DT-, and z-transform algorithms in this package scale logarithmically with the signal length $N$. As usual, the bond dimension keeps score.


## Authorship and maintenance

This repository supports the work reported in the (will update this after submitting on arxiv). 

- **Original author:** **Noufal Jaseem**
  Developed the original implementation and core workflow used to generate the results reported in the manuscript.

- **Lead maintainer / major contributor:** **Gauthameshwar S.**
Curated and structured the repository, improved code quality and usability, and implemented targeted efficiency optimizations. The optimized implementation was used to update the runtime figure in the manuscript.

- **Project supervision:** **Dario Poletti**
Provided scientific supervision and guidance for the project

## Contributing

Yes—please don’t hesitate to contribute and help develop this project further.
I am a poor PhD who learned how to create an open-source repository by starring popular packages, scrolling through Julia forums, and watching YouTube videos. I am very much not a full-fledged software developer.

If you find parts of this repo unpolished, unclear, or in need of improvement—and you think we could make it better together—please send a pull request and become a contributor 🫰🏻
Suggestions, refactors, documentation fixes, and performance ideas are all welcome.

## Citation

If you use QILaplace.jl in published work, please cite the repository and the paper that we will soon upload to arXiv.

---

Made with ❤️ using Julia, ITensors, and a healthy respect for bond-dimension growth.

---