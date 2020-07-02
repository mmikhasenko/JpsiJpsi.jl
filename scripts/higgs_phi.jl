using DrWatson
@quickactivate "JpsiJpsi"

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
using Plots
using DelimitedFiles
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
#
# sample creation
# S1_t = sample(500; H = Diagonal(fill(1/sqrt(3),3)))
# writedlm(datadir("simulations", "higgs_sample", "cosθ12ϕ.txt"),
#     [getproperty.(S1_t, :cosθ1) getproperty.(S1_t, :cosθ2) getproperty.(S1_t, :ϕ)]);

# sample loading
S1_t = let v = readdlm(datadir("simulations", "higgs_sample", "cosθ12ϕ.txt"))
    NamedTuple{(:cosθ1,:cosθ2,:ϕ)}.([v[i,:] for i in 1:size(v,1)])
end


using AlgebraPDF
using Statistics

Iϕ   = pdf((@. (ϕ;p)->1.0+p.β*cos(2ϕ)), p0=(β=0.3,), lims=(-π,π))
Icθi = pdf((@. (c;p)->1.0+p.ζ*(3*c^2-1)/2), p0=(ζ=0.3,), lims=(-1,1))
#
ϕv = getproperty.(S1_t, :ϕ)
frβ = fit_llh(ϕv, Iϕ; init_pars = [0.3])
#
cosθ1v = getproperty.(S1_t, :cosθ1)
cosθ2v = getproperty.(S1_t, :cosθ2)
frc = fit_llh(cosθ1v, I1; init_pars = [0.0])
frc = fit_llh(vcat(cosθ1v,cosθ2v), I1; init_pars = [0.0])
#

#                                                _|  _|
#    _|_|_|    _|_|_|  _|_|_|  _|_|    _|_|_|    _|      _|_|_|      _|_|_|
#  _|_|      _|    _|  _|    _|    _|  _|    _|  _|  _|  _|    _|  _|    _|
#      _|_|  _|    _|  _|    _|    _|  _|    _|  _|  _|  _|    _|  _|    _|
#  _|_|_|      _|_|_|  _|    _|    _|  _|_|_|    _|  _|  _|    _|    _|_|_|
#                                      _|                                _|
#                                      _|                            _|_|

@time βζs = [let Nev = 500
    S = sample(Nev; H = Diagonal(fill(1/sqrt(3),3)))
    ϕv = getproperty.(S, :ϕ)
    β = fit_llh(ϕv, Iϕ; init_pars = [0.3])[1]
    dv = getproperty.(S, :cosθ1)
    ζ = fit_llh(dv, Icθi; init_pars = [0.0])[1]
    (β=β, ζ=ζ)
end for _ in 1:500];

writedlm(datadir("simulations", "higgs_sample", "betazeta.txt"),
    [getproperty.(βζs, :β) getproperty.(βζs, :ζ)])

#            _|              _|      _|      _|
#  _|_|_|    _|    _|_|    _|_|_|_|_|_|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|
#  _|_|_|    _|    _|_|        _|_|    _|_|  _|  _|    _|    _|_|_|
#  _|                                                            _|
#  _|                                                        _|_|

fs = [(f=ϕ->Iϕ(ϕ,p=(β= 1/6,)), lab="0⁺ with b=d"),
      (f=ϕ->Iϕ(ϕ,p=(β=-1/4,)), lab="0⁻ with a=0"),
      (f=ϕ->Iϕ(ϕ,p=(β= 0.0,)), lab="1⁻"),
    ];

let Nbins=10
    Nev = length(S1_t)
    plot(size=(500,350), xlab="Δϕ", ylab="# events", xlim=(-π, π), leg=:top)
    [plot!(ϕ->f(ϕ) * (Nev/Nbins*(2π)), -π, π, lab=lab) for (f,lab) in fs]
    errorhist!(ϕv, bins=range(-π,π, length=Nbins+1), seriescolor=:black, lab="data")
    plot!(ylim=(0,70), leg = :bottomright)
    #
    plot!(inset=(1, bbox(0.15,0.5,0.3,0.35)))
    plot!(sp=2, xlab="β", ylab="# sample")
    stephist!(sp=2, getproperty.(βζs, :β), frame=:box, lab="", lc=:black)
    vline!(sp=2, [1/6], lab="", seriescolor=1, lw=2)
    vline!(sp=2, [  0], lab="", seriescolor=3, lw=2)
end
savefig(plotsdir("phi_higgs.pdf"))

βζs = let v = readdlm(datadir("simulations", "higgs_sample", "betazeta.txt"))
    NamedTuple{(:β,:ζ)}.([v[i,:] for i in 1:size(v,1)])
end

# let
#     plot(
#         stephist(getproperty.(βζs, :β), frame=:origin),
#         stephist(getproperty.(βζs, :ζ), frame=:origin)
#         )
#     # sqrt(cov(βs))
# end

#
# gs = [(f=cosθ->         (1)/2, lab="0+ with b=d"),
#       (f=cosθ->(1+1/2*cosθ)/2, lab="0-"),
#       (f=cosθ->(1-1/4*cosθ)/2, lab="1-"),
#     ];

gs = [(f=cosθ->Icθi(cosθ, p=(ζ=   0,)), lab="0+ with b=d"),
      (f=cosθ->Icθi(cosθ, p=(ζ= 1/2,)), lab="0-"),
      (f=cosθ->Icθi(cosθ, p=(ζ=-1/4,)), lab="1-"),
    ];

let Nbins=10
    Nev = length(S1_t)
    plot(size=(500,350), xlab="cosθ₁", ylab="# events", xlim=(-1, 1), leg=:top)
    [plot!(cosθ->f(cosθ) * (Nev/Nbins*(2)), -1, 1, lab=lab) for (f,lab) in gs]
    errorhist!(getproperty.(S1_t, :cosθ1), bins=range(-1,1, length=Nbins+1), seriescolor=:black, lab="data")
    plot!(ylim=(0,75), leg = :bottomright)
    #
    plot!(inset=(1, bbox(0.15,0.5,0.3,0.35)))
    plot!(sp=2, xlab="β", ylab="# sample")
    stephist!(sp=2, getproperty.(βζs, :ζ), frame=:box, lab="", lc=:black,
        bins=range(-0.3,0.55,length=20))
    vline!(sp=2, [   0], lab="", seriescolor=1, lw=2)
    vline!(sp=2, [ 1/2], lab="", seriescolor=2, lw=2)
    vline!(sp=2, [-1/4], lab="", seriescolor=3, lw=2)
end
savefig(plotsdir("costheta_higgs.pdf"))

#
let
    plot(size=(500,350),
        xlim=(-1.5,2.5), ylim=(-1.2,1.2),
        frame=:origin, leg=:bottomright)
    annotate!([(0.1,1.1,text("β", 10,:left)), (2.3,-0.05,text("ζ", 10, :top))])
    plot!(ξ->(2-ξ)/3, -1, 2, lab="group I", fill_between=0, seriescolor=1)
    plot!(ξ->(2-ξ)/3, -1, 2, lab="", lw=2, seriescolor=1)
    hline!([0.0], lab="", lc=:black)
    vline!([0.0], lab="", lc=:black)
    # plot!(ξ->-(1-2ξ)/3, -1, 1/2, lab="", fill_between=0, α=0.5, seriescolor=2)
    plot!(ξ->-(1-2ξ)/3, -1, 1/2, lab="group II", lw=8, seriescolor=2)
    scatter!([1/2], [0], lab="group III", ms=10, markerstrokewidth=0, markerstrokecolor=3, lc=:white, seriescolor=3)
    plot!(ξ->0, -1, 1/2, lab="group IV", lw=8, seriescolor=4)
    scatter!([1/2], [0], lab="", ms=10, markerstrokewidth=0, markerstrokecolor=3, seriescolor=3)
end

savefig(plotsdir("diagram_K.pdf"))

#
let
    plot(size=(500,350),
        xlim=(-1.3,0.8), ylim=(-0.3,0.3),
        yticks=([-1/4,0,1/4]),
        frame=:origin, leg=:bottomleft)
    annotate!([(0.05,0.27,text("β", 10,:left)), (0.7,-0.01,text("ζ", 10, :top))])
    plot!(ξ->(1+ξ)/6, -1, 1/2, lab="group I", fill_between=0, seriescolor=1)
    plot!(ξ->(1+ξ)/6, -1, 1/2, lab="", lw=2, seriescolor=1)
    hline!([0.0], lab="", lc=:black)
    vline!([0.0], lab="", lc=:black)
    # plot!(ξ->-(1/2+2ξ)/6, -1/4, 1/2, lab="", fill_between=0, α=0.5, seriescolor=2)
    plot!(ξ->-(1/2+2ξ)/6, -1/4, 1/2, lab="group II", lw=8, seriescolor=2)
    scatter!([-1/4], [0], lab="group III", ms=10, markerstrokewidth=0, markerstrokecolor=3, lw=0, seriescolor=3)
    plot!(ξ->0, -1/4, 1/2, lab="group IV", lw=8, seriescolor=4)
    scatter!([-1/4], [0], lab="", ms=10, markerstrokewidth=0, markerstrokecolor=3, lw=0, seriescolor=3)
end

savefig(plotsdir("diagram_mu.pdf"))

atan(-0.0,-0.0)
