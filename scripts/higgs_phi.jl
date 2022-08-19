using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
using Plots
using DelimitedFiles
using AlgebraPDF
using Statistics

theme(:wong)
include(joinpath("pyplot_settings.jl"))
#
@recipe function f(::Type{Val{:errorhist}}, x, y, z)
    h = Plots._make_hist((y,), plotattributes[:bins], normed=plotattributes[:normalize], weights=plotattributes[:weights])
    edge = collect(h.edges[1])
    centers = Plots._bin_centers(edge)
    dx = Plots.diff(h.edges[1]) / 2
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
    NamedTuple{(:cosθ1, :cosθ2, :ϕ)}.([v[i, :] for i in 1:size(v, 1)])
end


Iϕ = pdf((@. (ϕ; p) -> 1.0 + p.β * cos(2ϕ)), p0=(β=0.3,), lims=(-π, π))
Icθi = pdf((@. (c; p) -> 1.0 + p.ζ * (3 * c^2 - 1) / 2), p0=(ζ=0.3,), lims=(-1, 1))
#
ϕv = getproperty.(S1_t, :ϕ)
frβ = fit_llh(ϕv, Iϕ; init_pars=[0.3])[1]
#
cosθ1v = getproperty.(S1_t, :cosθ1)
cosθ2v = getproperty.(S1_t, :cosθ2)
frc = fit_llh(cosθ1v, Icθi; init_pars=[0.0])[1]
frc = fit_llh(vcat(cosθ1v, cosθ2v), I1; init_pars=[0.0])
#

#                                                _|  _|
#    _|_|_|    _|_|_|  _|_|_|  _|_|    _|_|_|    _|      _|_|_|      _|_|_|
#  _|_|      _|    _|  _|    _|    _|  _|    _|  _|  _|  _|    _|  _|    _|
#      _|_|  _|    _|  _|    _|    _|  _|    _|  _|  _|  _|    _|  _|    _|
#  _|_|_|      _|_|_|  _|    _|    _|  _|_|_|    _|  _|  _|    _|    _|_|_|
#                                      _|                                _|
#                                      _|                            _|_|

betazeta_filename = datadir("simulations", "higgs_sample", "betazeta.txt")
if !isfile(betazeta_filename)
    @time βζs = [
        let Nev = 500
            S = sample(Nev; H=Diagonal(fill(1 / sqrt(3), 3)))
            ϕv = getproperty.(S, :ϕ)
            β = fit_llh(ϕv, Iϕ; init_pars=[0.3])[1]
            dv = getproperty.(S, :cosθ1)
            ζ = fit_llh(dv, Icθi; init_pars=[0.0])[1]
            (β=β, ζ=ζ)
        end for _ in 1:500
    ]
    #
    writedlm(betazeta_filename,
        [getproperty.(βζs, :β) getproperty.(βζs, :ζ)])
end

# sqrt(cov(getproperty.(βζs, :β)))
# sqrt(cov(getproperty.(βζs, :ζ)))
βζs = let v = readdlm(betazeta_filename)
    NamedTuple{(:β, :ζ)}.([v[i, :] for i in 1:size(v, 1)])
end

#            _|              _|      _|      _|
#  _|_|_|    _|    _|_|    _|_|_|_|_|_|_|_|      _|_|_|      _|_|_|
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|
#  _|_|_|    _|    _|_|        _|_|    _|_|  _|  _|    _|    _|_|_|
#  _|                                                            _|
#  _|                                                        _|_|

fs = [(f=ϕ -> Iϕ(ϕ, p=(β=1 / 6,)), lab=L"0^+\,\,\mathrm{with}\,\,b=d"),
    (f=ϕ -> Iϕ(ϕ, p=(β=-1 / 4,)), lab=L"0^-"),
    (f=ϕ -> Iϕ(ϕ, p=(β=0.0,)), lab=L"1^-"),
];

let Nbins = 10
    Nev = length(S1_t)
    plot(size=(450, 300), xlab=L"\Delta \phi", ylab=L"\#\,\,\mathrm{events}", xlim=(-π, π), leg=:top)
    [plot!(ϕ -> f(ϕ) * (Nev / Nbins * (2π)), -π, π, lab=lab) for (f, lab) in fs]
    errorhist!(ϕv, bins=range(-π, π, length=Nbins + 1), seriescolor=:black, lab=L"\mathrm{data}")
    plot!(ylim=(0, 70), leg=:bottomright)
    #
    plot!(inset=(1, bbox(0.11, 0.51, 0.3, 0.35)))
    plot!(sp=2, xlab=L"\beta", ylab=L"\#\,\,\mathrm{samples}")
    stephist!(sp=2, getproperty.(βζs, :β), frame=:box, lab="", lc=:black,
        bins=range(-1 / 4, 0.32, length=35))
    vline!(sp=2, [1 / 6], lab="", seriescolor=1, lw=2)
    vline!(sp=2, [-1 / 4], lab="", seriescolor=2, lw=2)
    vline!(sp=2, [0], lab="", seriescolor=3, lw=2,
        xticks=([-1 / 4, 0, 1 / 6], ["-1/4", "0", "1/6"]))
end
savefig(plotsdir("phi_higgs.pdf"))

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

gs = [(f=cosθ -> Icθi(cosθ, p=(ζ=0,)), lab=L"0^{+}\,\,\mathrm{with}\,\,b=d"),
    (f=cosθ -> Icθi(cosθ, p=(ζ=1 / 2,)), lab=L"0^{-}"),
    (f=cosθ -> Icθi(cosθ, p=(ζ=-1 / 4,)), lab=L"1^{-}"),
];

let Nbins = 10
    Nev = length(S1_t)
    plot(size=(450, 300), xlab=L"\cos\theta_1", ylab=L"\#\,\,\mathrm{events}", xlim=(-1, 1), leg=:top)
    [plot!(cosθ -> f(cosθ) * (Nev / Nbins * (2)), -1, 1, lab=lab) for (f, lab) in gs]
    errorhist!(getproperty.(S1_t, :cosθ1), bins=range(-1, 1, length=Nbins + 1), seriescolor=:black, lab=L"\mathrm{data}")
    plot!(ylim=(0, 75), leg=:bottomright)
    #
    plot!(inset=(1, bbox(0.11, 0.51, 0.3, 0.35)))
    plot!(sp=2, xlab="", ann=(0.1, -14, text(L"\zeta", 10)), ylab=L"\#\,\,\mathrm{samples}")
    stephist!(sp=2, getproperty.(βζs, :ζ), frame=:box, lab="", lc=:black,
        bins=range(-1 / 4 - 0.04, 1 / 2 + 0.04, length=35), bottom_margin=0Plots.PlotMeasures.mm)
    vline!(sp=2, [0], lab="", seriescolor=1, lw=2)
    vline!(sp=2, [1 / 2], lab="", seriescolor=2, lw=2)
    vline!(sp=2, [-1 / 4], lab="", seriescolor=3, lw=2,
        xticks=([-1 / 4, 0, 1 / 2], ["-1/4", "0", "1/2"]))
end
savefig(plotsdir("costheta_higgs.pdf"))
