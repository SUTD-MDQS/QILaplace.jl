# rsvd.jl
# This file contains functions that performs rsvd based decomposition for a given ITensor
module RSVD

using ITensors
using LinearAlgebra
using Random

export rsvd

"""
    rsvd(A::ITensor, Linds...; k, p=10, q=0, random_seed, bondtag, verbose, cutoff, maxdim, mindim)

Low-allocation Randomized SVD (RSVD) for ITensor using native ITensor operations.
Returns (U, S, V) such that A ≈ U * S * V.

Optimizations:
  – Uses ITensor native `qr` and `svd` to avoid converting to dense Julia Arrays.
  – Minimizes allocations by keeping data in ITensor format.
  – Supports power iteration for better accuracy on decaying spectra.
"""
function rsvd(A::ITensor, Linds...; 
              k::Int=20, p::Int=10, q::Int=0, random_seed::Int=1234, 
              bondtag="Link,rsvd", verbose::Bool=false,
              cutoff::Float64=1e-15, maxdim::Int=k, mindim::Int=1)
    
    # 1. Setup indices
    Lis = commoninds(A, IndexSet(Linds...))
    Ris = uniqueinds(A, Lis)
    
    if length(Lis) == 0 || length(Ris) == 0
        error("In `rsvd`, left or right index set is empty. Left inds: $(Lis), right inds: $(Ris).")
    end

    # 2. Combiners (to treat A as a matrix)
    CL = combiner(Lis...)
    CR = combiner(Ris...)
    cL = combinedind(CL)
    cR = combinedind(CR)
    
    AC = A * CL * CR
    
    # 3. Random Projection
    # Rank target l = k + p (bounded by dimensions)
    l = min(k + p, dim(cL), dim(cR))
    
    Random.seed!(random_seed)
    α = Index(l, "rsvd_proj")
    Ω = random_itensor(eltype(A), cR, α)
    
    # Y = A * Ω
    Y = AC * Ω # (cL, α)
    
    # QR of Y -> Q
    # Q has indices (cL, q_link)
    Q, _ = qr(Y, cL; positive=true, tags="q_link")
    
    # 4. Power Iteration
    for iter in 1:q
        verbose && @info "RSVD power iteration $iter"
        # Z = A' * Q
        Z = dag(AC) * Q # (cR, q_link)
        Q_Z, _ = qr(Z, cR; positive=true, tags="q_link")
        
        # Y = A * Q_Z
        Y = AC * Q_Z # (cL, q_link)
        Q, _ = qr(Y, cL; positive=true, tags="q_link")
    end
    
    # 5. Form small matrix B = Q' * A
    B = dag(Q) * AC 
    
    q_link = commonind(B, Q) # The link index connecting Q and B
    
    # 6. SVD of small matrix B. We also pass truncation parameters here.
    U_small, S, V = svd(B, q_link; 
                        cutoff=cutoff, maxdim=maxdim, mindim=mindim,
                        lefttags=bondtag, righttags=bondtag)
    
    # 7. Form full U = Q * U_small
    U = Q * U_small # (cL, u_link)
    
    # 8. Uncombine
    U = U * dag(CL)
    V = V * dag(CR)
    
    return U, S, V
end

end # module RSVD
