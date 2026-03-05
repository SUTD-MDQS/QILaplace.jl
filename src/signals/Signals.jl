# src/signals/Signals.jl
module Signals

using ITensors, Random, Printf
using LinearAlgebra: normalize, normalize!

export generate_signal

######################################## BASIC SIGNALS ########################################

# --- multi-sine --------------------------------------------------------
function _generate_signal(
    ::Val{:sin},
    n::Int,
    dt::Float64,
    freq::Float64;
    phase::Float64=0.0,
    noise_level::Float64=0.0,
    kwargs...,
)
    jvals = 0:(2^n-1)
    return [sin(freq * dt * j + phase) + noise_level * randn() for j in jvals]
end

function _generate_signal(
    ::Val{:sin},
    n::Int,
    dt::Float64,
    freq::Vector{Float64};
    phase::Vector{Float64}=zeros(length(freq)),
    noise_level::Float64=0.0,
    kwargs...,
)
    jvals = 0:(2^n-1)
    length(freq) == length(phase) ||
        throw(ArgumentError("Frequency and phase vectors must be of the same length."))
    return [
        sum(sin(ω * dt * j + φ) for (ω, φ) in zip(freq, phase)) + noise_level * randn() for
        j in jvals
    ]
end

# --- pure random ------------------------------------------------------------
function _generate_signal(::Val{:random}, n::Int, dt::Float64; seed::Int=1234, kwargs...)
    rng = MersenneTwister(seed)
    return randn(rng, 2^n)
end

# --- multi-sin with exponential envelope -----------------------------------
function _generate_signal(
    ::Val{:sin_decay},
    n::Int,
    dt::Float64,
    freq::Float64;
    decay_rate::Float64,
    phase::Float64=0.0,
    kwargs...,
)
    jvals = 0:(2^n-1)
    return [sin(freq * dt * j + phase) * exp(-decay_rate * dt * j) for j in jvals]
end

function _generate_signal(
    ::Val{:sin_decay},
    n::Int,
    dt::Float64,
    freq::Vector{Float64};
    decay_rate::Vector{Float64},
    phase::Union{Vector{Float64},Nothing}=nothing,
    kwargs...,
)
    jvals = 0:(2^n-1)
    length(freq) == length(decay_rate) ||
        throw(ArgumentError("Frequency and decay_rate vectors must be of the same length."))

    if phase === nothing
        return [
            sum(sin(ω * dt * j) * exp(-λ * dt * j) for (ω, λ) in zip(freq, decay_rate)) for
            j in jvals
        ]
    else
        length(freq) == length(phase) ||
            throw(ArgumentError("Frequency and phase vectors must be of the same length."))
        return [
            sum(
                sin(ω * dt * j + φ) * exp(-λ * dt * j) for
                (ω, λ, φ) in zip(freq, decay_rate, phase)
            ) for j in jvals
        ]
    end
end

##################################### Root dispatch function #########################################
"""
    generate_signal(n; kind=:sin, dt=nothing, freq=nothing, kwargs...)

Generate a length-`2^n` real signal of a specified type.

# Arguments
- `n::Int`: The exponent determining the length of the signal. The signal will have `2^n` points.

# Keyword Arguments
- `kind::Symbol`: The type of signal to generate. Defaults to `:sin`. Supported kinds include:
    - `:sin`: A sinusoidal signal.
    - `:cos`: A cosine signal.
    - `:decay`: An exponentially decaying sinusoidal signal.
    - `:noise`: A purely random noise signal.
- `dt::Union{Nothing, Float64}`: The time step between samples. If `nothing`, it is automatically 
   computed from `freq` as `dt = 1 / (freq * 2^n)`. Defaults to `nothing`.
- `freq::Union{Nothing, Float64, Vector{Float64}}`: The frequency or frequencies of the signal. 
   Can be a scalar or a vector of frequencies. Defaults to `2π` if not provided.

# kwargs (kind-dependent)
- `phase::Float64`: The phase offset of the signal in radians. Applicable to `:sin` and `:cos` kinds. 
   Defaults to `0.0`.
- `decay_rate::Float64`: The rate of exponential decay. Only applicable to the `:decay` kind. 
   A higher value results in faster decay. Defaults to `1.0`.
- `noise_level::Float64`: The amplitude of random noise added to the signal. Applicable to all kinds. 
   Defaults to `0.0` (no noise).
- `seed::Int`: The random seed for reproducibility of noise generation. Applicable when `noise_level > 0`. 
   Defaults to `nothing` (random seed).

# Returns
- `signal::Vector{Float64}`: A real-valued vector of length `2^n` representing the generated signal.

# Examples
```julia
n = 5

# simple sine wave with frequency 1.0
x = generate_signal(n, kind=:sin, freq=1.0) 

# multi-frequency decaying sine wave
x = generate_signal(n, kind=:sin_decay, freq=[1.0, 2.0], decay_rate=[0.1, 0.2])
```
"""
function generate_signal(
    n::Int;
    kind::Symbol=:sin,
    dt::Union{Float64,Nothing}=nothing,
    freq::Union{Real,Vector{<:Real},Nothing}=nothing,
    kwargs...,
)

    # 1. Handle Random (no freq needed)
    if kind == :random
        dt_val = isnothing(dt) ? 1.0 : Float64(dt)
        return _generate_signal(Val(:random), n, dt_val; kwargs...)
    end

    # 2. Handle Frequency-based signals
    freq_val = isnothing(freq) ? 2π : freq
    freq_typed = freq_val isa Real ? Float64(freq_val) : Float64.(freq_val)

    # Default dt logic (preserving your formulas)
    if isnothing(dt)
        if freq_typed isa Vector
            f_max = maximum(abs, freq_typed)
            dt_val = f_max == 0 ? 1.0 : (2π / f_max / n)
        else
            dt_val = freq_typed == 0 ? 1.0 : (2 * freq_typed / n)
        end
    else
        dt_val = Float64(dt)
    end

    if kind == :sin
        return _generate_signal(Val(:sin), n, dt_val, freq_typed; kwargs...)
    elseif kind == :sin_decay
        return _generate_signal(Val(:sin_decay), n, dt_val, freq_typed; kwargs...)
    else
        throw(
            ArgumentError(
                "Unsupported signal kind: $kind. Supported kinds are :sin, :sin_decay, :random.",
            ),
        )
    end
end

end # module Signals
