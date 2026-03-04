# # Discrete Fourier Transform Tutorial
#
# This walkthrough shows how to use QILaplace.jl to compute the DFT.

using QILaplace, ITensors

# ## Setting up the signal
# Create an 8-point signal:

n = 8
x = generate_signal(n, kind=:sin, freq=1.0)

# ## Constructing the 'SignalMPS' from the data vector

sites = [Index(2, "site-$i") for i in 1:n]
psi_test, _ = signal_mps(x)

# ## Constructing the QFT circuit
qft_mpo = build_qft_mpo(psi_test; cutoff=1e-14, maxdim=1000)

# Ensure your MPS and MPO have the same site indices before applying. (note!)
# ## Doing the Fourier Transform
result = apply(qft_mpo, psi_test)
