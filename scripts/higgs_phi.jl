using DrWatson
@quickactivate "JpsiJpsi"

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
using Plots
theme(:wong)
#
@recipe function f(::Type{Val{:errorhist}}, x, y, z)
    h = Plots._make_hist((y,), plotattributes[:bins], normed = plotattributes[:normalize], weights = plotattributes[:weights])
    edge = collect(h.edges[1])
    centers = Plots._bin_centers(edge)
    dx = Plots.diff(h.edges[1])/2
    x := centers
    y := h.weights
    seriestype := :scatter
    @series begin
        seriestype := :xerror
        xerror --> dx
        yerror --> sqrt.(h.weights)
        (centers, h.weights)
    end
    @series begin
        seriestype := :yerror
        yerror --> sqrt.(h.weights)
        (centers, h.weights)
    end
    ()
end
@shorthands errorhist
#

fs = [(f=ϕ->(1+cos(2*ϕ)/6)/2π, lab="0⁺ with b=d"),
      (f=ϕ->(1-cos(2*ϕ)/4)/2π, lab="0⁻ with a=0"),
      (f=ϕ->1/2π, lab="1⁻"),
    ];
S1_t = sample(500; H = Diagonal(fill(1/sqrt(3),3)))

let Nbins=10
    Nev = length(S1_t)
    plot(size=(500,350), xlab="Δϕ", ylab="# events", xlim=(-π, π), leg=:top)
    [plot!(ϕ->f(ϕ) * (Nev/Nbins*(2π)), -π, π, lab=lab) for (f,lab) in fs]
    errorhist!(getproperty.(S1_t, :ϕ), bins=range(-π,π, length=Nbins+1), seriescolor=:black, lab="data")
    plot!(ylim=(0,70), leg = :bottomright)
end
savefig(plotsdir("phi_higgs.pdf"))


gs = [(f=cosθ->(1)/2, lab="0+ with b=d"),
      (f=cosθ->(1+1/2*cosθ)/2, lab="0-"),
      (f=cosθ->(1-1/4*cosθ)/2, lab="1-"),
    ];

let Nbins=10
    Nev = length(S1_t)
    plot(size=(500,350), xlab="cosθ₁", ylab="# events", xlim=(-1, 1), leg=:top)
    [plot!(cosθ->f(cosθ) * (Nev/Nbins*(2)), -1, 1, lab=lab) for (f,lab) in gs]
    errorhist!(getproperty.(S1_t, :cosθ1), bins=range(-1,1, length=Nbins+1), seriescolor=:black, lab="data")
    plot!(ylim=(0,75), leg = :bottomright)
end
savefig(plotsdir("costheta_higgs.pdf"))

#
let
    plot(size=(500,350),
        xlim=(-1.5,2.5), ylim=(-1.2,1.2),
        frame=:origin, leg=:bottomright)
    annotate!([(0.1,1.1,text("β", 10,:left)), (2.3,-0.05,text("ξ", 10, :top))])
    plot!(ξ->(2-ξ)/3, -1, 2, lab="group I", fill_between=0, α=0.5, seriescolor=1)
    plot!(ξ->(2-ξ)/3, -1, 2, lab="", lw=2, seriescolor=1)
    # plot!(ξ->-(1-2ξ)/3, -1, 1/2, lab="", fill_between=0, α=0.5, seriescolor=2)
    plot!(ξ->-(1-2ξ)/3, -1, 1/2, lab="group II", lw=6, seriescolor=2)
    plot!(ξ->0, -1, 1/2, lab="group IV", lw=6, seriescolor=4)
    scatter!([1/2], [0], lab="group III", ms=10, markerstrokewidth=0, l=nothing, seriescolor=3)
end

savefig(plotsdir("diagram_K.pdf"))

#
let
    plot(size=(500,350),
        xlim=(-1.3,0.8), ylim=(-0.3,0.3),
        yticks=([-1/4,0,1/4]),
        frame=:origin, leg=:bottomleft)
    annotate!([(0.05,0.27,text("β", 10,:left)), (0.7,-0.01,text("ξ", 10, :top))])
    plot!(ξ->(1+ξ)/6, -1, 1/2, lab="group I", fill_between=0, α=0.5, seriescolor=1)
    plot!(ξ->(1+ξ)/6, -1, 1/2, lab="", lw=2, seriescolor=1)
    # plot!(ξ->-(1/2+2ξ)/6, -1/4, 1/2, lab="", fill_between=0, α=0.5, seriescolor=2)
    plot!(ξ->-(1/2+2ξ)/6, -1/4, 1/2, lab="group II", lw=6, seriescolor=2)
    plot!(ξ->0, -1/4, 1/2, lab="group IV", lw=6, seriescolor=4)
    scatter!([-1/4], [0], lab="group III", ms=10, markerstrokewidth=0, l=nothing, seriescolor=3)
end

savefig(plotsdir("diagram_mu.pdf"))
