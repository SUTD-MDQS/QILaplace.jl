```@meta
CurrentModule = QILaplace
```

# QILaplace

QILaplace.jl is a Julia library for quantum-inspired signal transforms on classical hardware. It combines Matrix Product State (MPS) data representations with Matrix Product Operator (MPO) circuit compression to evaluate Fourier- and Laplace-type transforms at scales where dense methods become impractical. Instead of materializing full vectors or matrices, QILaplace works directly with compressed tensor networks, so the workflow stays memory-efficient and controllable through truncation tolerance and bond dimension. 

The package currently provides three transform families in one consistent API: a Quantum Fourier Transform, a real-axis Laplace (damping) transform, and a full complex discrete Laplace (z-transform). Around these, it includes signal conversion to MPS using SVD and randomized SVD support for scalable preprocessing.

## Installation

You can install QILaplace by running the following command in the Julia REPL:
```julia
using Pkg;Pkg.add("QILaplace.jl");
```
If you are more of an interactive coder, you can run the same by going to the pkg mode: 
```julia
] add QILaplace
```

This package is built on ITensors.jl, which requires Julia 1.10 or above for a stable performance. So it is recommended to have Julia 1.10+

## Start Here

- [Core Concepts](core_concepts.md): learn the tensor-network ideas behind the package.
- [Signal Tutorial](tutorials/signal.md): start with signal encoding and compression.
- [DFT Tutorial](tutorials/dft.md): see the QFT-style Fourier workflow.
- [DT Tutorial](tutorials/dt.md): understand the damping transform.
- [zT Tutorial](tutorials/zt.md): follow the full discrete Laplace / z-transform workflow.
- [Benchmarking](benchmarking.md): inspect runtime and scaling behavior.
- [API Reference](api.md): browse the public interface.

## Further Learning

If you are interested in the concept/theory of laplace, and fourier transforms, feel free to check out these videos by 3B1B. 

```@raw html
<div style="display:flex; gap:24px; flex-wrap:wrap; margin:20px 0;">
  <div style="max-width:320px;">
    <a href="https://www.youtube.com/watch?v=spUNpyF58BY" target="_blank" rel="noopener">
      <img src="https://img.youtube.com/vi/spUNpyF58BY/hqdefault.jpg" alt="3Blue1Brown Fourier transform video" style="width:100%; border-radius:12px; border:1px solid #ccc;">
    </a>
    <p><strong>3Blue1Brown:</strong> A visual introduction to Fourier transforms.</p>
  </div>
  <div style="max-width:320px;">
    <a href="https://www.youtube.com/watch?v=j0wJBEZdwLs" target="_blank" rel="noopener">
      <img src="https://img.youtube.com/vi/j0wJBEZdwLs/hqdefault.jpg" alt="3Blue1Brown Laplace transform video" style="width:100%; border-radius:12px; border:1px solid #ccc;">
    </a>
    <p><strong>3Blue1Brown:</strong> Visual intuition for Laplace-transform ideas.</p>
  </div>
</div>

```@autodocs
Modules = [QILaplace]
```

## Further Resources
- [TensorNetwork.org](https://tensornetwork.org/) for a solid introduction to MPS, MPOs, and tensor-network methods.
