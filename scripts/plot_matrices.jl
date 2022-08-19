using Plots
using JpsiJpsi
using PartialWaveFunctions
using QuadGK

# plotting
let J = 2
    plot(heatmap.(nHs[J+1], c=:balance, clim=(-1, 1))...,
        layout=grid(1, length(nHs[J+1])), size=(200 * length(nHs[J+1]), 190), colorbar=false,
        xaxis=false, yaxis=false)
end
let J = 1
    plot(heatmap.(nHs[J+1], c=:balance, clim=(-1, 1))...,
        layout=grid(1, length(nHs[J+1])), size=(200 * length(nHs[J+1]), 190), colorbar=false,
        xaxis=false, yaxis=false)
end
let J = 0
    plot(heatmap.(nHs[J+1], c=:balance, clim=(-1, 1))...,
        layout=grid(1, length(nHs[J+1])), size=(200 * length(nHs[J+1]), 190), colorbar=false,
        xaxis=false, yaxis=false)
end

function Hλλ(λ1, λ2; LS, J, j1=1, j2=1)
    L, S = LS
    Δλ = λ1 - λ2
    (isodd(j2 - λ2) ? -1 : 1) *
    CG(j1, λ1, j2, -λ2, S, Δλ) *
    CG(L, 0, S, Δλ, J, Δλ)
end


let
    e(m, m′) =
        quadgk(z -> 3 * sum(wignerd(1, m, ξ1, z) * wignerd(1, m′, ξ1, z) for ξ1 in [-1, 1]), 0, 1)[1] / 2
    [e(m, m′) for m in -1:1, m′ in -1:1]
end


# f(m,m′) = quadgk(z->3*sum(wignerd(1,m,ξ1,z)*wignerd(1,m′,ξ1,z) for ξ1 in [-1,1]), -1, 1)[1]/2
f(m, m′) =
    quadgk(z -> 3 * sum(wignerd(1, m, ξ1, z) * wignerd(1, m′, ξ1, z) +
                        wignerd(1, m, ξ1, -z) * wignerd(1, m′, ξ1, -z) for ξ1 in [-1, 1]), 0, 1)[1] / 2

F = [f(m, m′) for m in -1:1, m′ in -1:1];

g(m, m′) = quadgk(z -> 3 * wignerd(1, m, 0, z) * wignerd(1, m′, 0, z), -1, 1)[1] / 2

[g(m, m′) for m in -1:1, m′ in -1:1]


h(m, m′) =
    quadgk(z -> 3 * sum(wignerd(1, m, ξ1, z) * wignerd(1, m′, ξ1, z) -
                        wignerd(1, m, ξ1, -z) * wignerd(1, m′, ξ1, -z) for ξ1 in [-1, 1]), 0, 1)[1] / 2
#
S = [h(m, m′) for m in -1:1, m′ in -1:1] .* sqrt(2);

k(m, m′) =
    quadgk(z -> 3 * sum(wignerd(1, m, ξ1, z) * wignerd(1, m′, ξ1, z) +
                        wignerd(1, m, ξ1, 1 - z) * wignerd(1, m′, ξ1, 1 - z) for ξ1 in [-1, 1]), 0, 1)[1] / 2
#
K = [k(m, m′) for m in -1:1, m′ in -1:1];
K


S
I(H, F1, F2) =
    sum(H[λ1, λ2] * F1[λ1, λ1′] * H[λ1′, λ2′] * F2[λ2, λ2′] * kronecker(λ1 - λ2, λ1′ - λ2′)
        for λ1 in 1:3, λ2 in 1:3, λ1′ in 1:3, λ2′ in 1:3)
#
