# rsvd.jl
# This file contains functions that performs rsvd based decomposition for a given ITensor
module RSVD

using ITensors
using LinearAlgebra
using Random

export rsvd

"""
    rsvd(A::ITensor, Linds...; k=20, p=10, q=0, random_seed=1234, bondtag="Link,rsvd",
         verbose=false, cutoff=1e-15, maxdim=k, mindim=1) -> (U, S, V)

Randomised SVD for an `ITensor`, returning factors `A ≈ U * S * V` where `S`
is diagonal. Uses native ITensor `qr`/`svd` internally to avoid dense matrix
conversions.

# Arguments
- `A`       — the ITensor to decompose.
- `Linds`   — indices that should end up on the left factor `U`.

# Keyword Arguments
- `k::Int=20`          — target rank (number of kept singular values).
- `p::Int=10`          — oversampling; the projection uses `k + p` random vectors.
- `q::Int=0`           — number of power iterations (improves accuracy on slowly decaying spectra).
- `random_seed::Int`   — seed for the random projection.
- `cutoff::Float64=1e-15` — singular value truncation threshold.
- `maxdim::Int=k`      — hard cap on bond dimension.
- `mindim::Int=1`      — minimum bond dimension to keep.
- `bondtag`            — tag string applied to the SVD bond index.
- `verbose::Bool`      — print power-iteration progress.

# When to use
Use `:rsvd` (via `signal_mps(x; method=:rsvd, k=...)`) when the signal is very large
and you want a fast, low-rank approximation instead of a full SVD.
"""
function rsvd(
    A::ITensor,
    Linds...;
    k::Int=20,
    p::Int=10,
    q::Int=0,
    random_seed::Int=1234,
    bondtag="Link,rsvd",
    verbose::Bool=false,
    cutoff::Float64=1e-15,
    maxdim::Int=k,
    mindim::Int=1,
)

    # 1. Setup indices
    Lis = commoninds(A, IndexSet(Linds...))
    Ris = uniqueinds(A, Lis)

    if length(Lis) == 0 || length(Ris) == 0
        error(
            "In `rsvd`, left or right index set is empty. Left inds: $(Lis), right inds: $(Ris).",
        )
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
    U_small, S, V = svd(
        B,
        q_link;
        cutoff=cutoff,
        maxdim=maxdim,
        mindim=mindim,
        lefttags=bondtag,
        righttags=bondtag,
    )

    # 7. Form full U = Q * U_small
    U = Q * U_small # (cL, u_link)

    # 8. Uncombine
    U = U * dag(CL)
    V = V * dag(CR)

    return U, S, V
end

end # module RSVD
