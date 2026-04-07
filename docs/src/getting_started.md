# Core Concepts

To effectively use QILaplace.jl, it is helpful to understand how classical signals are mapped onto tensor network structures. This library leverages Matrix Product States (MPS) and Matrix Product Operators (MPO) to perform transformations—like the Laplace and Fourier transform—with logarithmic scaling.

## Matrix Product States (MPS) and Operators (MPO)

At its core, QILaplace.jl does not treat a signal as a long, flat array of numbers. Instead, it reshapes that data into a high-dimensional tensor and compresses it into an MPS.

### What are MPS and MPO?
A Matrix Product State (MPS) is a mathematical factorization that breaks down a massive tensor into a chain of smaller, interconnected local tensors. While originally popularized in quantum physics to represent many-body wavefunctions, in signal processing, an MPS serves as a highly efficient compressed data structure.

- Efficiency: A signal of size $N = 2^n$ usually requires $O(N)$ memory. An MPS can represent the same signal using only $O(n \cdot \chi^2)$ parameters, where $\chi$ is the bond dimension.
- Correlation: The bond dimension $\chi$ represents the "information complexity" or the amount of correlation between different segments of the signal. For many physical and mathematical signals, $\chi$ remains small even as the signal size $N$ grows exponentially.

A Matrix Product Operator (MPO) is the operator equivalent of an MPS. It allows us to apply linear transformations—such as the Fourier or Laplace transform—directly to the compressed MPS without ever needing to decompress the data back into a flat array.

```@raw html
<div style="background-color: #eef9f0; border: 1px solid #75d689ff; border-radius: 10px; padding: 14px 16px; margin: 20px 0;">
    <div style="color: #000000; font-weight: 800; font-size: 1.02rem; letter-spacing: 0.03em; margin-bottom: 6px;;">💡 Further Reading</div>
    <div style="color: #1f4a2a;">For a deeper dive into the physics-style origins and mathematical theory of MPS, we recommend exploring <a href="https://tensornetwork.org/mps/" target="_blank" rel="noopener">TensorNetwork.org</a> or the seminal review by <a href="https://arxiv.org/abs/1008.3477" target="_blank" rel="noopener">Schollwöck (2011)</a>.</div>
</div>
```

## Compressing Data: The Quantics Representation

The bridge between a standard 1D signal and a compressed MPS is a process called Binary Encoding (often referred to in literature as the Quantics representation).

### Binary-encoding of Signal
To process a signal of length $N = 2^n$, we treat the index of each data point as a binary string of length $n$. We then reshape the vector into an $n$-dimensional tensor of shape $(2, 2, \dots, 2)$.
In QILaplace.jl, we follow a big-endian encoding convention. For a signal vector $x$ where we want to access the element $x_j$:

1. We represent the index $j$ as a bit-string: $j = (b_1, b_2, \dots, b_n)_2$.
2. The first bit $b_1$ corresponds to the "slowest" changing index (leftmost in the tensor chain).
3. The last bit $b_n$ corresponds to the "fastest" changing index (rightmost).

In the animation below, you can see an array of 8 entries ($x_0$ to $x_7$) getting rearranged into a 3-dimensional cube. Each axis of this cube corresponds to one of the three indices of our new tensor. By looping through the tensor indices (using the $1, 2$ convention), we can map every point in the high-dimensional space back to its original location in the flat signal array.

```@raw html
<div style="display: flex; justify-content: center;">
    <video controls width="800" style="border: 1px solid #3B4252; border-radius: 15px; margin: 20px 0;">
        <source src="../animations/assets/TensorifyVector.mp4" type="video/mp4">
        Your browser does not support the video tag.
    </video>
</div>
```

## MPS Conversion Algorithms

Transforming the "tensorized" signal into a compressed MPS is the most computationally intensive phase of the entire pipeline. This stage serves as the primary bottleneck because it requires decomposing the exponential information of the full signal into a chain of local correlations. Regardless of the specific algorithm used, the conversion process relies on the decay of Singular Values. As we decompose the signal, we encounter a spectrum of singular values at each bond connecting the tensor sites.

- Compression: We "truncate" the representation by keeping only the most significant singular values.
- The Threshold ($\tau$): Users can define a relative cutoff $\tau$. Any singular value smaller than this threshold is discarded.
- Accuracy-Compression Trade-off: A smaller $\tau$ leads to higher fidelity but larger bond dimensions ($\chi$), while a larger $\tau$ achieves massive compression at the cost of some numerical precision.

To perform this decomposition, QILaplace.jl provides two primary algorithms, each meticulously optimized for different hardware constraints and signal complexities.

### 1. The Standard SVD (Sequential Sweep)
The standard approach involves a sequential sweep from one end of the tensor chain to the other.

- Mechanism: At each site, the algorithm performs a Singular Value Decomposition (SVD) to split the current tensor into a local MPS site and a remainder that is passed to the next site.
- Cost: The bottleneck occurs at the central bond of the chain (where the matricized tensor is largest). For an $n$-qubit system, the complexity at the center is $O(2^{3n/2})$, as it requires computing the full singular value spectrum of a $2^{n/2} \times 2^{n/2}$ matrix.
- Best for: Small to medium $n$ where the exact spectrum is needed for high-fidelity representation.

```@raw html
<div style="display: flex; justify-content: center;">
    <video controls width="800" style="border: 1px solid #3B4252; border-radius: 15px; margin: 20px 0;">
        <source src="../animations/assets/SVDRMPSConversion.mp4" type="video/mp4">
        Your browser does not support the video tag.
    </video>
</div>
```

### 2. Randomized SVD (Divide and Conquer)
The Randomized SVD (RSVD) algorithm uses a "divide-and-conquer" strategy to significantly speed up the conversion.

- Mechanism: Instead of sweeping linearly, RSVD divides the tensor at the middle into a "left" and "right" block. It then "conquers" by iteratively splitting these blocks into single-site tensors.
- Approximation: Unlike standard SVD, RSVD finds an approximation of the top-$k$ singular values using random projections. This avoids the need to find the full spectrum.
- Cost: By targeting only the relevant $k$ singular values (where $k \approx \chi$), the complexity at the central bond is reduced to $O(k \cdot 2^{n})$, offering a massive speedup for large signals.
- Best for: Large-scale signals (e.g., $n > 20$) where the signal is known to have a low rank ($k \ll N$).

```@raw html
<div style="display: flex; justify-content: center;">
    <video controls width="800" style="border: 1px solid #3B4252; border-radius: 15px; margin: 20px 0;">
        <source src="../animations/assets/RSVDRMPSConversion.mp4" type="video/mp4">
        Your browser does not support the video tag.
    </video>
</div>
```

This dual-algorithm approach allows QILaplace.jl to balance the trade-off between absolute numerical precision and the ability to process exponentially large data sets on standard hardware.

## Quantum-Inspired Signal Processing

The "Quantum-Inspired" nature of this library is defined by how we construct our operators. In traditional quantum computing, a circuit is a sequence of unitary gates acting on qubits. In QILaplace.jl, we reinterpret these circuits as Matrix Product Operators (MPOs) that perform linear transformations directly on an MPS.

### From Circuits to Compressed Operators
By "zipping" circuit gates into an MPO, we can compress an entire transformation into a compact form with low bond dimension.

- Unitary Transforms ($Q\hat{F}T$): The Quantum Fourier Transform is the foundation of this approach. As introduced by [Chen, Stoudenmire, and White (2023)](https://arxiv.org/abs/2210.08468), the QFT circuit can be highly compressed into an MPO where the bond dimension does not increase with the number of qubits $n$. This allows spectral analysis on exponentially large data in logarithmic time.

- Non-Unitary Transforms ($\hat{DT}$): QILaplace.jl extends this by incorporating non-unitary maps. The Discrete Laplace Transform requires exponential damping, which we implement via a Damping Transform ($\hat{DT}$). This circuit uses non-unitary damping gates that are likewise compressed into an efficient MPO.

By operating on classical hardware, we gain a unique "digital advantage." We can merge the $\hat{DT}$ and $Q\hat{F}T$ into a single, combined $z\hat{T}$ MPO. This unified operator allows us to probe the complex $z$-plane and identify the poles and zeros of a signal at scales reaching $M = 2^{60}$ points—far exceeding the limits of standard FFT-based methods.

```@raw html
<div style="background-color: #fff8e8; border: 1px solid #e7a747ff; border-radius: 10px; padding: 14px 16px; margin: 20px 0;">
    <div style="color: #000000; font-weight: 800; font-size: 1.02rem; letter-spacing: 0.03em; margin-bottom: 6px;">💡 Why Quantum-Inspired?</div>
    <div style="color: #6a4a1b;">We borrow the "algorithmic structure" of quantum gates but implement them as compressed classical tensors. This allows us to execute non-unitary maps that are often difficult for real quantum computers to handle, while maintaining the exponential scaling benefits of quantum algorithms.</div>
</div>
```

### Quantum Fourier Transform Circuit

The Quantum Fourier Transform (QFT) is one of the many foundational quantum algorithms that demonstrate an exponential speedup. While the classical FFT requires $O(N \log N)$ operations for a signal of size $N$, the QFT algorithm theoretically operates in $O(\log^2 N)$ time. In our context, we don't execute gates on a physical device; instead, we represent the entire circuit structure as a single, static MPO.

To turn the QFT circuit into an efficient MPO, we use the Zip-Up algorithm followed by a sweep of the Orthogonality Center (OC) with truncation, as introduced by [Chen, Stoudenmire, and White (2023)](https://arxiv.org/abs/2210.08468). Instead of multiplying all gates together (which would lead to an exponentially large matrix), the algorithm "zips" the gates into a chain of local tensors and moves the OC cleverly while truncating so that the changes are global.

The process follows a specific sequence as seen in the animation below:

**Gate Combination:** Tensors representing individual gates are combined sequentially across the qubits. These gates lie in the same qubit line.

**QR & SVD Truncation:** At each step, a Singular Value Decomposition (SVD) is performed on the bonds. This identifies the "entanglement" the operator introduces between sites and zips up two adjacent tensor trains.

**Compression:** We move down the chain, truncating small singular values at the OC.

```@raw html
<div style="display: flex; justify-content: center;">
    <video controls width="800" style="border: 1px solid #3B4252; border-radius: 15px; margin: 20px 0;">
        <source src="../animations/assets/CircuitCompression.mp4" type="video/mp4">
        Your browser does not support the video tag.
    </video>
</div>
```

Theoretical Guarantee: The efficiency of this compression relies on the fact that the singular values $\sigma_k$ for the QFT and Laplace operators decay exponentially with distance. Specifically, the singular values at a bond often follow the relationship:


$$\sigma_k \sim e^{-\alpha k}$$


where $k$ is the index of the singular value and $\alpha > 0$ is a decay constant. This ensure that our QFT circuit is represented in an MPO that does not scale with system size. Without this guarantee, we could attempt to compress any arbitrary quantum circuit, but the bond dimensions would grow uncontrollably with qubit size and circuit depth.
