```@meta
CurrentModule = QILaplace
```

# QILaplace
(An introductory statement about QILaplace.jl here)

(Some detailed description of what this does)

(Features this package has: DFT, inverse DFT, real exponent Laplace transform/damping transform, Discrete Laplace transform/z transform)

(When QILaplace can do good signal processing and when does it not)

Links to the other sections of this page. Go to the Introduction to know the basic concepts used in this repository, tutorials to know more about blablabla, and benchmarking to see how well the quantum-inspired algorithms work on your computer!

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
