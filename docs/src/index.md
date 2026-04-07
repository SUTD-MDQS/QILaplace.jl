```@meta
CurrentModule = QILaplace
```

# QILaplace
(An introductory statement about QILaplace.jl here)

(Some detailed description of what this does)

(Features this package has: DFT, inverse DFT, real exponent Laplace transform/damping transform, Discrete Laplace transform/z transform)

(When QILaplace can do good signal processing and when does it not)

Links to the other sections of this page. Explore [Core Concepts](core_concepts.md) to understand the theoretical foundations of the library, head to the [Tutorials](tutorials/signal.md) for hands-on coding examples, and check [Benchmarking](benchmarking.md) to see how these quantum-inspired algorithms perform on your hardware.

## Tutorials

- [Signal Encoding and Compression Tutorial](tutorials/signal.md)
- [Discrete Fourier Transform Tutorial](tutorials/dft.md)


```@autodocs
Modules = [QILaplace]
```

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
