using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using FileIO
using Plots
using LaTeXStrings
using Statistics

theme(:wong, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    xlim=(:auto, :auto), ylim=(:auto, :auto))
# include(joinpath("pyplot_settings.jl"))

function get_fit_results(filename)
    f = JLD2.jldopen(filename)
    JLD2.read(f, "H1_t"), JLD2.read(f, "TS")
end

wl = load(datadir("simulations", "groups_cross_testing", "Natt=5_Nev=500_Nsampl=500_gen_group=1_higgs_1m11.jld2"));
fr_0p = wl["fit_results"];

let bins = range(-1, -0.7, length=100)
    plot(size=(450, 300), leg=:topleft,
        xlab=L"\mathcal{L}(G\{\hat{h}\})",
        ylab=L"\#\,\,\mathrm{sample\,\,entries}")
    stephist!(-map(v -> v[1].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab="", seriescolor=1, α=0.4, fill_between=0)
    stephist!(-map(v -> v[1].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab=L"\mathrm{fit\,\,by\,\,group-I}", seriescolor=1, lw=1.0)
    stephist!(-map(v -> v[2].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab="", seriescolor=2, α=0.4, fill_between=0)
    stephist!(-map(v -> v[2].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab=L"\mathrm{fit\,\,by\,\,group-II}", seriescolor=2, lw=1.0)
    stephist!(-map(v -> v[3].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab="", seriescolor=3, α=0.4, fill_between=0)
    stephist!(-map(v -> v[3].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab=L"\mathrm{fit\,\,by\,\,group-III}", seriescolor=3, lw=1.0)
    stephist!(-map(v -> v[4].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab="", seriescolor=4, α=0.4, fill_between=0)
    stephist!(-map(v -> v[4].LLH, fr_0p) ./ wl["Nev"], bins=bins, lab=L"\mathrm{fit\,\,by\,\,group-IV}", seriescolor=4, lw=1.0)
end
savefig(plotsdir("llh_testing_higgs.pdf"))

_, fr_1p = get_fit_results(datadir("simulations", "groups_cross_testing", "results_J3_Nev500_Natt5_Nsamples501.jld2"));
const Nev_1m = 500;

let bins = range(-0.84, -0.75, length=100)
    plot(size=(500, 350))
    stephist!(-map(v -> v[1].LLH, fr_1p) ./ Nev_1m, bins=bins, lab="group-I")
    stephist!(-map(v -> v[2].LLH, fr_1p) ./ Nev_1m, bins=bins, lab="group-II")
    stephist!(-map(v -> v[3].LLH, fr_1p) ./ Nev_1m, bins=bins, lab="group-III")
    stephist!(-map(v -> v[4].LLH, fr_1p) ./ Nev_1m, bins=bins, lab="group-IV")
end

let bins = range(-30 / 500, 90 / 500, length=50)
    plot(size=(450, 300), xlab=L"\mathrm{TS}(0^+/1^-)", ylab=L"\#\,\,\mathrm{sample\,\,entries}")
    stephist!(-map(v -> v[1].LLH - v[3].LLH, fr_0p) ./ Nev_1m, bins=bins, lab="", seriescolor=1, fill_between=0, alpha=0.4)
    stephist!(-map(v -> v[1].LLH - v[3].LLH, fr_0p) ./ Nev_1m, bins=bins, lab=L"0^+\,\,\mathrm{generated}", seriescolor=1, lw=1.0)
    stephist!(-map(v -> v[1].LLH - v[3].LLH, fr_1p) ./ Nev_1m, bins=bins, lab="", seriescolor=3, fill_between=0, alpha=0.4)
    stephist!(-map(v -> v[1].LLH - v[3].LLH, fr_1p) ./ Nev_1m, bins=bins, lab=L"1^-\,\,\mathrm{generated}", seriescolor=3, lw=1.0)
end
savefig(plotsdir("TS_0p_vs_1m.pdf"))

m1, σ1 = let v = -map(v -> v[1].LLH - v[3].LLH, fr_0p) ./ Nev_1m
    mean(v), sqrt(cov(v))
end

m2, σ2 = let v = -map(v -> v[1].LLH - v[3].LLH, fr_1p) ./ Nev_1m
    mean(v), sqrt(cov(v))
end

m1, σ1
m2, σ2
(m1 - m2) / sqrt(σ1^2 + σ2^2)
