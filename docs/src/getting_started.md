# Introduction

## Introduction to Matrix Product States, and Matrix Product Operators
This is a section giving an overview of MPS, MPO and why they are efficient in compressing certain types of data

## Compressing data into MPS
This section will give a summary of the binary encoding, how this type of encoding is used in quantum computers to store classical bits in qubits. The advantage is that we require only log(N) qubits to store an array of size N. We can then convert this data tensor into an MPS which compresses the data by only storing the relevant correlation information in its bonds. 

## Quantum-Inspired Signal Processing
This is a section telling about how the early works of signal processing involved compressing the QFT circuit efficiently in the form of MPO with bond dimensions that do not scale with system size. Then a motivation on how non-unitary gates can also be incorporated in this framework and can be represented in the MPO form, thus leading to a truly quantum-inspired algorithm that has an advantage in its digital form. 

