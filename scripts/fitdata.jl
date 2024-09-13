using Plots
using DelimitedFiles
using Statistics
using Optim
using AlgebraPDF
using LinearAlgebra
using RecipesBase

Base.length(f::FunctionWithParameters) = 1
Base.iterate(f::FunctionWithParameters) = (f, nothing)
Base.iterate(f::FunctionWithParameters, state) = nothing

theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    xlim=(:auto, :auto), ylim=(:auto, :auto))
# 
cms_file = joinpath("data", "dijpsi_cms.dlm")
data = readdlm(cms_file)

struct diJpsiData{T}
    data::T
end

@recipe function plot(x::diJpsiData)
    binwidth = mean(diff(data[:, 1]))
    xerror --> binwidth .* one.(x.data[:, 1])
    yerror --> sqrt.(x.data[:, 2])
    seriestype --> :scatter
    markersize --> 3
    markerstyle --> :d
    xguide --> "m [GeV]"
    yguide --> "dσ/dm"
    (x.data[:, 1], x.data[:, 2])
end


plot(diJpsiData(data), label="CMS data")


# 
const mth = 6.18  # GeV
continuum =
    FunctionWithParameters(
        (m; p) -> exp(-p.α * (m - mth)) / p.α * (1 + p.a * (m - mth)),
        Ext(α=0.85, a=0.0))
# 
BW1 = FunctionWithParameters(
    (m; p) -> (p.m1 * p.Γ1)^2 / abs2(p.m1^2 - m^2 - 1im * p.m1 * p.Γ1),
    Ext(m1=6.5, Γ1=0.3))
# 
BW2 = FunctionWithParameters(
    (m; p) -> (p.m2 * p.Γ2)^2 / abs2(p.m2^2 - m^2 - 1im * p.m2 * p.Γ2),
    Ext(m2=6.9, Γ2=0.1))
# 
BW3 = FunctionWithParameters(
    (m; p) -> (p.m3 * p.Γ3)^2 / abs2(p.m3^2 - m^2 - 1im * p.m3 * p.Γ3),
    Ext(m3=7.3, Γ3=0.05))
# 

ph = FunctionWithParameters((m; p) -> sqrt(m - mth), ∅)
SB = FSum([continuum, BW1, BW2, BW3] .* ph, Ext(fb=170.0, fs1=100.0, fs2=50.0, fs3=20.0))


plot!(SB, 6.3, 9)
χ2 = ChiSq(SB, data[:, 1], data[:, 2], sqrt.(data[:, 2])) |>
     m -> fixpars(m, (:m1, :Γ3))

fit_model = let
    loss = Base.Fix1(χ2, ())
    p0 = pars(χ2)
    δp0 = (α=0.009, a=1.0,
        m1=0.1, Γ1=0.01, m2=0.1, Γ2=0.01, m3=0.1, Γ3=0.01,
        fb=1.0, fs1=1.0, fs2=1.0, fs3=1.0)
    # 
    invH_stepsize = Diagonal(p2v(δp0, χ2)) .+ eps()
    initial_invH = x -> invH_stepsize
    res = optimize(loss, p2v(p0, χ2), BFGS(; initial_invH))
    fit_pars = v2p(Optim.minimizer(res), χ2)
    updatepars(χ2.f, fit_pars)
end

let
    plot()
    plot!(fit_model, mth, 9, label="Fit model", l=(2, :red))
    plot!(fit_model[1], mth, 9, label="Continuum")
    plot!(fit_model[2], mth, 9, fill=0, α=0.4, label="Resonance 1")
    plot!(fit_model[3], mth, 9, fill=0, α=0.4, label="Resonance 2")
    plot!(fit_model[4], mth, 9, fill=0, α=0.4, label="Resonance 3")
    plot!(diJpsiData(data), label="CMS data", c=:black)
end
