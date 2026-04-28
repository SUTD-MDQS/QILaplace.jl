# API

This page documents the public interface of QILaplace.jl.

## MPS Representation

```@docs
QILaplace.SignalMPS
QILaplace.ZTMPS
QILaplace.Mps.PairCore
```

### Construction
```@docs
QILaplace.generate_signal
QILaplace.signal_mps
QILaplace.signal_ztmps
```

### Extraction & Indexing
```@docs
QILaplace.coefficient
QILaplace.mps_to_vector
QILaplace.siteindices
QILaplace.bondindices
```

### Manipulation
```@docs
QILaplace.canonicalize!
QILaplace.compress!
QILaplace.update_site!
QILaplace.update_bond!
QILaplace.Mps.norm
```

## MPO Representation

```@docs
QILaplace.SingleSiteMPO
QILaplace.PairedSiteMPO
```

## Transformers

```@docs
QILaplace.build_qft_mpo
QILaplace.build_dt_mpo
QILaplace.build_zt_mpo
```

## Linear Algebra Utilities

```@docs
QILaplace.apply
QILaplace.RSVD.rsvd
```
