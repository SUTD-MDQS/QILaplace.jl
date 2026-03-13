# Introduction

## Introduction to Matrix Product States, and Matrix Product Operators
This is a section giving an overview of MPS, MPO and why they are efficient in compressing certain types of data

## Compressing data into MPS
This section will give a summary of the binary encoding, how this type of encoding is used in quantum computers to store classical bits in qubits. The advantage is that we require only log(N) qubits to store an array of size N. We can then convert this data tensor into an MPS which compresses the data by only storing the relevant correlation information in its bonds. 


### Binary-encoding of Signal
```@raw html
<video controls width="900" style="border: 1px solid #3B4252; border-radius: 15px;">
	<source src="../animations/assets/TensorifyVector.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>
```

### SVD to MPS Conversion
The SVD (Singular Value Decomposition) to MPS conversion is a key algorithmic procedure for transforming arbitrary tensors into the canonical MPS form. This process involves iteratively splitting tensors using SVD, where each singular value matrix is absorbed into the neighboring tensor. The animation below shows how a 4-index tensor is progressively decomposed through three SVD operations, resulting in a chain of left-orthogonal matrices connected by bond indices χ. The bond dimensions (χ) represent the entanglement complexity at each bond and can be truncated to control the accuracy-compression trade-off.

```@raw html
<video controls width="900" style="border: 1px solid #3B4252; border-radius: 15px;">
	<source src="../animations/assets/SVDRMPSConversion.mp4" type="video/mp4">
	Your browser does not support the video tag.
</video>
```

## Quantum-Inspired Signal Processing
This is a section telling about how the early works of signal processing involved compressing the QFT circuit efficiently in the form of MPO with bond dimensions that do not scale with system size. Then a motivation on how non-unitary gates can also be incorporated in this framework and can be represented in the MPO form, thus leading to a truly quantum-inspired algorithm that has an advantage in its digital form. 

