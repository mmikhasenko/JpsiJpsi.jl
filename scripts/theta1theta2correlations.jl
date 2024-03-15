using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using JpsiJpsi
using Parameters
using Plots
using LaTeXStrings
using PartialWaveFunctions

theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    xlim=(:auto, :auto), ylim=(:auto, :auto))


function I4μ_intϕ(vars; H=error("helicity coupling"))
    @unpack cosθ1, cosθ2 = vars
    return real(sum(
        wignerd(1, λ1, ξ1, cosθ1) * wignerd(1, λ2, ξ2, cosθ2) *
        wignerd(1, λ1, ξ1, cosθ1) * wignerd(1, λ2, ξ2, cosθ2) *
        abs2(H[λ1+2, λ2+2])
        for ξ1 in [-1, 1], ξ2 in [-1, 1],
        λ1 in -1:1, λ2 in -1:1))
end


function LS2λλ(λ1, λ2; LS, J, j1=1, j2=1)
    L, S = LS
    Δλ = λ1 - λ2
    (isodd(j2 - λ2) ? -1 : 1) *
    clebschgordan(j1, λ1, j2, -λ2, S, Δλ) *
    clebschgordan(L, 0, S, Δλ, J, Δλ)
end

I4μ_intϕ(vars, J::Int, LS::Tuple{Int,Int}) =
    I4μ_intϕ(vars; H=[LS2λλ(λ1, λ2; LS=LS, J) for λ1 in -1:1, λ2 in -1:1])

let (L, S) = (0, 0), J = 0
    cv = range(-1, 1, length=100)
    v = [I4μ_intϕ((cosθ1=c1, cosθ2=c2), J, (L, S))
         for (c1, c2) in Iterators.product(cv, cv)]
    v ./= sum(v) / (100 / 2)^2
    plot(size=(450, 400),
        xlab=L"\cos\,\theta_1", ylab=L"\cos\,\theta_2", title="J = $J: (L=$L, S=$S)")
    heatmap!(cv, cv, v, c=:viridis)
end



# LS couplings
const LS_for_J0 = [(0, 0), (2, 2)];
const LS_for_J1 = [(0, 1), (2, 1), (2, 2)];
const LS_for_J2 = [(0, 2), (2, 0), (2, 1), (2, 2), (4, 2)];

const JHs = [[[LS2λλ(λ1, λ2; J=i - 1, LS) for λ1 = -1:1, λ2 = -1:1] for LS in LSs]
             for (i, LSs) in enumerate([LS_for_J0, LS_for_J1, LS_for_J2])]
const nJHs = [[h ./ sqrt(intensity(h)) for h in Hs] for Hs in JHs]

let
    pv = []
    for (J, LSs) in (0 => LS_for_J0, 1 => LS_for_J1, 2 => LS_for_J2)
        for (L, S) in LSs
            cv = range(-1, 1, length=100)
            v = [I4μ_intϕ((cosθ1=c1, cosθ2=c2), J, (L, S))
                 for (c1, c2) in Iterators.product(cv, cv)]
            v ./= sum(v) / (100 / 2)^2
            plot(size=(450, 400),
                xlab=L"\cos\,\theta_1", ylab=L"\cos\,\theta_2", title="J = $J: (L=$L, S=$S)")
            heatmap!(cv, cv, v, c=:viridis, colorbar=false)
            push!(pv, plot!())
        end
    end
    plot(pv..., size=(300 * 4, 280 * 3))
end
savefig(joinpath("plots", "map_JLS.pdf"))
