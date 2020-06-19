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

fs = [(f=ϕ->(1+cos(2*ϕ)/6)/2π, lab="\"0⁺\""),
      (f=ϕ->(1-cos(2*ϕ)/6)/2π, lab="\"0⁻\""),
      (f=ϕ->1/2π, lab="1⁺/1⁻"),
    ];
S1_t = sample(1000; H = Diagonal(fill(1/sqrt(3),3)))

let Nbins=15
    Nev = length(S1_t)
    plot(size=(500,350), xlab="Δϕ", ylab="# events", xlim=(-π, π), leg=:top)
    [plot!(ϕ->f(ϕ) * (Nev/Nbins*(2π)), -π, π, lab=lab) for (f,lab) in fs]
    errorhist!(getproperty.(S1_t, :ϕ), bins=range(-π,π, length=Nbins+1), seriescolor=:black, lab="data")
end
savefig(plotsdir("phi_higgs.pdf"))
