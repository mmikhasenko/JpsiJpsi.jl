using DrWatson
@quickactivate "JpsiJpsi"

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
#

fs = [(f=ϕ->(1+cos(2*ϕ)/6)/2π, lab="\"0⁺\""),
      (f=ϕ->(1-cos(2*ϕ)/6)/2π, lab="\"0⁻\""),
      (f=ϕ->1/2π, lab="1⁺/1⁻"),
    ];
S1_t = sample(1000; H = Diagonal(fill(1/sqrt(3),3)))

let
    Nev = length(S1_t)
    plot(size=(500,350), xlab="ϕ", ylab="# events", leg=:bottomright, xlim=(-π, π))
    [plot!(ϕ->f(ϕ) * (Nev/10*(2π)), -π, π, lab=lab) for (f,lab) in fs]
    stephist!(getproperty.(S1_t, :ϕ), bins=range(-π,π, length=11), seriescolor=:black, lab="data")
end
savefig(plotsdir("phi_higgs.pdf"))
