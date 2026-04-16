# # Discrete Fourier Transform Tutorial
#
# This walkthrough shows how to use `QILaplace.jl` to compute the DFT.
#
# The Quantum Fourier Transform (QFT) is the unitary transform underlying several
# quantum algorithms. In modern form it was formalized in early quantum computing
# work (notably Coppersmith, 1994) and popularized by Shor's factoring algorithm.
# In this tutorial, we use the same linear map on a *classical* tensor-network state.
#
# For a register of $n$ qubits ($N = 2^n$ basis states), the QFT is defined by
#
# ```math
# \mathrm{QFT}_N\,|x\rangle = \frac{1}{\sqrt{N}}\sum_{k=0}^{N-1}
# e^{2\pi ixk/N}|k\rangle,
# ```
#
# where $x,k\in\{0,\dots,N-1\}$ and $|x\rangle$ denotes the computational basis
# state corresponding to the binary encoding of index $x$: $|x\rangle = |x_1x_2\dots x_n\rangle$.
#
# In `QILaplace.jl`, we first encode a length-$N$ signal into an MPS,
# $\sum_{x=0}^{N-1} a_x |x\rangle$, then apply a compressed QFT MPO built from
# Hadamard and controlled-phase gates.
#
# ![QFT circuit for $n=4$ qubits](../animations/assets/qft_circuit.png)
#
# *Figure 1. Four-qubit QFT circuit in product form: Each wire receives a Hadamard gate followed by controlled-phase gates ($P_{ij}$) from less significant wires. The unswapped output is in bit-reversed order.*
#
# In the product-form circuit, outputs come out in **bit-reversed order** relative
# to standard DFT indexing: if $k = k_{n-1}\dots k_1k_0$, the unswapped output is
# indexed by $\operatorname{rev}(k)=k_0k_1\dots k_{n-1}$.
#
# This is why we either (i) append explicit SWAP gates at the end of the QFT
# circuit, (ii) sample the MPS by keeping the bit-reversed ordering in mind, or
# (iii) account for reversal when converting the transformed MPS to a dense
# vector for comparison with FFT-based reference implementations.
#
# In this QFT circuit, the Hadamard Gate ($H$) and the Phase gate ($P_{ij}$) are defined as:
# ```math
# H = \frac{1}{\sqrt{2}} \begin{pmatrix} 1 & 1 \\ 1 & -1 \end{pmatrix} \quad
# P_{ij} = \begin{pmatrix} 1 & 0 \\ 0 & e^{2\pi i / 2^{j-i+1}} \end{pmatrix}
# ```
# Controlled-phase action:
# - If the control qubit is $|0\rangle$, nothing happens to the target.
# - If the control qubit is $|1\rangle$, the target picks up the phase defined by $P_{ij}$.
#

using QILaplace, ITensors
using FFTW, LinearAlgebra

# ## Setting up the signal
# Create a $2^n$-point signal (here $n=4$, so $N=16$):

n = 4
x = generate_signal(n, kind=:sin, freq=1.0)

# ## Constructing the SignalMPS

sites = [Index(2, "site-$i") for i in 1:n]
psi_test, x_norm = signal_mps(x)

#
# Coefficient-level validation and signal-compression diagnostics are now covered
# in the dedicated `signal.jl` tutorial so that this page stays focused on DFT/QFT.

# ## Constructing the QFT circuit
qft_mpo = build_qft_mpo(psi_test; cutoff=1e-14, maxdim=100)

# Ensure MPS and MPO use exactly the same site indices before applying.
qft_mpo.sites == psi_test.sites # true #hide

# ## Performing the Fourier Transform
psi_qn = apply(qft_mpo, psi_test)

#
# You can also simply multiply the QFT MPO with the SignalMPS using `qft_mpo * psi_test`.
#
# To compare the results of our transform, we use `FFTW.jl` as the reference to verify our transform results. 
# FFTW conventions are:
#
# ```math
# \mathrm{fft}(x)_k = \sum_{x=0}^{N-1} x_x\,e^{-2\pi i xk/N},
# \qquad
# \mathrm{bfft}(x)_k = \sum_{x=0}^{N-1} x_x\,e^{+2\pi i xk/N}.
# ```
#
# Our QFT convention uses the $+2\pi i$ phase and includes $1/\sqrt{N}$, so for
# normalized input $\hat{x}=x/\|x\|_2$ we expect
#
# ```math
# \mathrm{QFT}_N\hat{x} = \frac{\mathrm{bfft}(\hat{x})}{\sqrt{N}}.
# ```

# We now need a helper that contracts a `SignalMPS` back to a dense vector. The
# option `rev=true` applies output-index bit reversal, so the returned vector is
# in the usual DFT ordering.

function bitreverse_int(v::Integer, n::Int) #hide
	bits_lsb = digits(v; base=2, pad=n) #hide
	return sum(bits_lsb[i] << (n - i) for i in 1:n) #hide
end #hide
#
function mps_to_vector(psi::SignalMPS; rev::Bool=false) #hide
	nsites = length(psi.sites) #hide
	N = 2^nsites #hide

	T = prod(psi.data) #hide
	arr = Array(T, reverse(psi.sites)...) #hide
	vec_qn = reshape(arr, N) #hide

	if !rev #hide
		return vec_qn #hide
	end #hide

	perm = [bitreverse_int(k, nsites) + 1 for k in 0:(N - 1)] #hide
	return vec_qn[perm] #hide
end #hide
nothing #hide

# `mps_to_vector` is practical for small to moderate $n$ and debugging. For very
# large systems, dense reconstruction is exponentially expensive in memory/time.
# We recommend sampling the spectrum directly from the MPS form in such cases. 
#
# Next we compare QILaplace output with FFTW and sample a few indices.

N = length(x)
x_hat = x / x_norm

qft_qn = mps_to_vector(psi_qn; rev=false)
qft_fn = mps_to_vector(psi_qn; rev=true)
fftw_ref = bfft(x_hat) / sqrt(N)

println("\nQFT coefficients in circuit order (unswapped):")
println(round.(qft_qn; digits=5))

println("\nQFT coefficients in FFT order (bit-reversed):")
println(round.(qft_fn; digits=5))

comparison_error = norm(qft_fn - fftw_ref)
@show comparison_error

# Explicit sampling for k = 1 and k = 7.
# n=4:
# - k=1  -> binary 0001 -> reversed 1000 -> index 8 in 0-based -> 9 in Julia.
# - k=7  -> binary 0111 -> reversed 1110 -> index 14 in 0-based -> 15 in Julia.
k1 = 1
idx_fn_k1 = k1 + 1
idx_qn_k1 = 9
qft_val_k1 = round(qft_qn[idx_qn_k1]; digits=5)
fftw_val_k1 = round(fftw_ref[idx_fn_k1]; digits=5)
@show k1 qft_val_k1 fftw_val_k1 abs(qft_val_k1 - fftw_val_k1)

k2 = 7
idx_fn_k2 = k2 + 1
idx_qn_k2 = 15
qft_val_k2 = round(qft_qn[idx_qn_k2]; digits=5)
fftw_val_k2 = round(fftw_ref[idx_fn_k2]; digits=5)
@show k2 qft_val_k2 fftw_val_k2 abs(qft_val_k2 - fftw_val_k2)

# For a moderately larger signal, we compare the QFT spectrum against FFTW and plot the
# absolute error on a secondary (right) y-axis with a distinct linestyle.

using Plots, LaTeXStrings

n_big = 8
x_big = generate_signal(n_big, kind=:sin, freq=[4.0, 17.0], phase=[0.0, 0.3])
psi_big, x_big_norm = signal_mps(x_big)

qft_big_mpo = build_qft_mpo(psi_big; cutoff=1e-12, maxdim=1000)
psi_big_qn = apply(qft_big_mpo, psi_big)

qft_big = mps_to_vector(psi_big_qn; rev=true)
fftw_big = bfft(x_big / x_big_norm) / sqrt(length(x_big))
abs_err_big = abs.(qft_big .- fftw_big)

N_big = length(x_big)

# For this generated signal, we can predict where peaks should appear, and what the DC value $(at \omega=0$) is analytically.
#
# The signal generator uses
#
# ```math
# x_j = \sum_r \sin(\Omega_r j + \phi_r),\quad \Omega_r = \omega_r\,dt,
# ```
#
# and for vector frequencies it sets
#
# ```math
# dt = \frac{2\pi}{\omega_{\max}\,n}.
# ```
#
# With `freq=[4,17]`, `phase=[0,0.3]`, and `n=8`, we get
#
# ```math
# \Omega_1 = \frac{\pi}{17},\qquad \Omega_2 = \frac{\pi}{4}.
# ```
#
# So the spectrum should show symmetric peaks near
# $\omega\approx\pm\pi/17$ and $\omega=\pm\pi/4$.
#
# The DC value (at $\omega=0$) for the normalized transform
# $\mathrm{bfft}(x/\|x\|_2)/\sqrt{N}$ equals
#
# ```math
# X(0) = \frac{1}{\sqrt{N}\,\|x\|_2}\sum_{j=0}^{N-1}x_j.
# ```
#
# For each sinusoid, the finite sum is
#
# ```math
# S(\Omega,\phi)=\sum_{j=0}^{N-1}\sin(\Omega j+\phi)
# =\frac{\sin\!\left(\frac{N\Omega}{2}\right)
# \sin\!\left(\phi+\frac{(N-1)\Omega}{2}\right)}{\sin\!\left(\frac{\Omega}{2}\right)}.
# ```

freq_big = [4.0, 17.0]
phase_big = [0.0, 0.3]
dt_big = (2π / maximum(abs, freq_big) / n_big)
Ω_big = freq_big .* dt_big

@show Ω_big Ω_big ./ π

function sine_sum_closed_form(Ω::Real, ϕ::Real, N::Int)
	return sin(N * Ω / 2) * sin(ϕ + (N - 1) * Ω / 2) / sin(Ω / 2)
end

dc_pred = sum(sine_sum_closed_form(Ω, ϕ, N_big) for (Ω, ϕ) in zip(Ω_big, phase_big)) / (x_big_norm * sqrt(N_big))

# Shift spectra so frequency runs from approximately -$\pi$ to $\pi$.
qft_big_shift = fftshift(qft_big)
fftw_big_shift = fftshift(fftw_big)

ω_axis = (2π / N_big) .* collect(-div(N_big, 2):(div(N_big, 2) - 1))
peak_marks = [-π / 4, -π / 17, π / 17, π / 4];
#
xtick_vals = sort(unique(vcat([-π, -π / 2], peak_marks, [0.0, π / 2, π]))) #hide
xtick_labels = map(xtick_vals) do ω #hide
	if isapprox(ω, -π / 4; atol=1e-10) #hide
		L"-\frac{\pi}{4}" #hide
	elseif isapprox(ω, -π / 17; atol=1e-10) #hide
		L"-\frac{\pi}{17}" #hide
	elseif isapprox(ω, π / 17; atol=1e-10) #hide
		L"\frac{\pi}{17}" #hide
	elseif isapprox(ω, π / 4; atol=1e-10) #hide
		L"\frac{\pi}{4}" #hide
	elseif isapprox(ω, -π / 2; atol=1e-10) #hide
		L"-\frac{\pi}{2}" #hide
	elseif isapprox(ω, π / 2; atol=1e-10) #hide
		L"\frac{\pi}{2}" #hide
	elseif isapprox(ω, -π; atol=1e-10) #hide
		L"-\pi" #hide
	elseif isapprox(ω, π; atol=1e-10) #hide
		L"\pi" #hide
	else #hide
		"" #hide
	end #hide
end #hide
nothing #hide
p = plot( #hide
	ω_axis, #hide
	abs.(qft_big_shift); #hide
    label="|QILaplace QFT|", #hide
    linewidth=2, #hide
	xlabel=L"\omega", #hide
    ylabel="Magnitude", #hide
	xticks=(xtick_vals, xtick_labels), #hide
	legend=:topleft, #hide
); #hide
plot!(p, ω_axis, abs.(fftw_big_shift); label="|FFTW bfft|", linewidth=2, linestyle=:dash); #hide
vline!(p, peak_marks; color=:grey, linestyle=:dashdot, linewidth=1.0, label=false); #hide
#
p_err = twinx(p); #hide
abs_err_big_shift = fftshift(abs_err_big); #hide
err_max = maximum(abs_err_big_shift); #hide
err_ylim = max(2 * err_max, eps()); #hide
plot!( #hide
    p_err, #hide
	ω_axis, #hide
	abs_err_big_shift; #hide
    label="|error|", #hide
	linewidth=1.5, #hide
    linestyle=:dot, #hide
    color=:grey, #hide
    ylabel="Absolute error", #hide
	legend=:topright, #hide
); #hide
ylims!(p_err, (0.0, err_ylim)); #hide
#
plot_path = joinpath(@__DIR__, "..", "assets", "dft_spectrum_comparison.svg"); #hide
mkpath(dirname(plot_path)); #hide
savefig(p, plot_path); #hide
nothing #hide

zero_idx = findfirst(==(0.0), ω_axis)
dc_numeric = isnothing(zero_idx) ? NaN + NaN * im : fftw_big_shift[zero_idx]
@show dc_pred dc_numeric abs(dc_pred - dc_numeric)

# The generated spectrum comparison plot is embedded below.
#
# ```@raw html
# <img src="../../assets/dft_spectrum_comparison.svg" alt="QFT spectrum vs FFTW with absolute error">
# ```
#
# *Figure 2. Shifted spectrum comparison in angular frequency $\omega\in[-\pi,\pi)$ for $n=8$: the QILaplace QFT and FFTW reference overlap at the expected peak locations near $\omega\approx\pm\pi/17$ and $\omega=\pm\pi/4$, while the dotted error curve (right axis) remains small throughout the band.*
# 
# Currently we don't have the inverse QFT available in `QILaplace.jl`, but it is straightforward to implement since QFT is a unitary, hence invertible transform. If you want this support as well, feel free to leave a request in our [main GutHub repo](https://github.com/SUTD-MDQS/QILaplace.jl/issues) :)
