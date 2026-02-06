# test/test_aqua.jl
# This file tests whether QILaplace.jl passes the standards of Aqua.jl
using Aqua, QILaplace

@testset "Aqua.jl" begin
    Aqua.test_all(QILaplace)
end
