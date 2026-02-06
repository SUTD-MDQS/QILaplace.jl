using QILaplace.Signals

@testset "Signals.jl: basic functionality" begin
    n = 10

    # 1. Basic Sin
    s1 = generate_signal(n, kind=:sin, freq=1.0)
    @test length(s1) == 2^n
    @test s1 isa Vector{Float64}

    # 2. Vector Freq Sin
    s2 = generate_signal(n, kind=:sin, freq=[1.0, 2.0])
    @test length(s2) == 2^n

    # 3. Random
    s3 = generate_signal(n, kind=:random)
    @test length(s3) == 2^n

    # 4. Sin Decay (scalar)
    s4 = generate_signal(n, kind=:sin_decay, freq=1.0, decay_rate=0.1)
    @test length(s4) == 2^n

    # 5. Sin Decay (vector)
    s5 = generate_signal(n, kind=:sin_decay, freq=[1.0, 2.0], decay_rate=[0.1, 0.2])
    @test length(s5) == 2^n

    # 6. Error handling
    @test_throws ArgumentError generate_signal(
        n, kind=:sin_decay, freq=[1.0], decay_rate=[0.1, 0.2]
    ) # Mismatch
    @test_throws ArgumentError generate_signal(n, kind=:unsupported)
end

@testset "Signals.jl: analytical correctness and randomness" begin
    # simple analytic sine (no noise)
    n = 2
    dt = 1.0
    freq = 1.0
    phase = 0.3
    s = generate_signal(n, kind=:sin, dt=dt, freq=freq, phase=phase, noise_level=0.0)
    expected = [sin(freq * dt * j + phase) for j in 0:(2 ^ n - 1)]
    @test all(isapprox.(s, expected; atol=1e-12, rtol=0))

    # multi-frequency with phases
    s2 = generate_signal(
        n, kind=:sin, dt=dt, freq=[1.0, 2.0], phase=[0.0, 0.5], noise_level=0.0
    )
    expected2 = [
        sum(sin(ω * dt * j + φ) for (ω, φ) in zip([1.0, 2.0], [0.0, 0.5])) for
        j in 0:(2 ^ n - 1)
    ]
    @test all(isapprox.(s2, expected2; atol=1e-12, rtol=0))

    # sin_decay scalar
    s3 = generate_signal(n, kind=:sin_decay, dt=dt, freq=2.0, decay_rate=0.1, phase=0.0)
    expected3 = [sin(2.0 * dt * j) * exp(-0.1 * dt * j) for j in 0:(2 ^ n - 1)]
    @test all(isapprox.(s3, expected3; atol=1e-12, rtol=0))

    # sin_decay vector and default phase
    s4 = generate_signal(n, kind=:sin_decay, dt=dt, freq=[1.0, 2.0], decay_rate=[0.1, 0.2])
    expected4 = [
        sum(sin(ω * dt * j) * exp(-λ * dt * j) for (ω, λ) in zip([1.0, 2.0], [0.1, 0.2]))
        for j in 0:(2 ^ n - 1)
    ]
    @test all(isapprox.(s4, expected4; atol=1e-12, rtol=0))

    # random seed deterministic
    nrand = 4
    r1 = generate_signal(nrand, kind=:random, seed=42)
    r2 = generate_signal(nrand, kind=:random, seed=42)
    @test r1 == r2

    # integer freq should work like float
    sint = generate_signal(n, kind=:sin, dt=dt, freq=1)
    sfloat = generate_signal(n, kind=:sin, dt=dt, freq=1.0)
    @test all(isapprox.(sint, sfloat; atol=1e-12, rtol=0))
end
