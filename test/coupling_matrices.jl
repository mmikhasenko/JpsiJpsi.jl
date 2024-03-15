using JpsiJpsi
using Test


let
    # J = 0
    [Hλλ(λ1, λ2; LS=(0, 0), J=0) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(2, 2), J=0) for λ1 in -1:1, λ2 in -1:1]
    # J = 1
    [Hλλ(λ1, λ2; LS=(0, 1), J=1) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(2, 1), J=1) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(2, 2), J=1) for λ1 in -1:1, λ2 in -1:1]
    # J = 2
    [Hλλ(λ1, λ2; LS=(0, 2), J=2) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(2, 0), J=2) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(2, 1), J=2) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(2, 2), J=2) for λ1 in -1:1, λ2 in -1:1]
    [Hλλ(λ1, λ2; LS=(4, 2), J=2) for λ1 in -1:1, λ2 in -1:1]
    #
    @test true
end
#
@test A4μ((cosθ1=0.1, ϕ1=0.2, cosθ2=0.3, ϕ2=0.4); LS=(0, 0), J=0) != 0

@testset "Orthogonality" begin
    @test prod([hi × hj for hi in ngHs[1], hj in ngHs[3]] .== 0)
    @test prod([hi × hj for hi in ngHs[1], hj in ngHs[2]] .== 0)
    @test prod([hi × hj for hi in ngHs[3], hj in ngHs[2]] .== 0)
end