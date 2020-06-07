using DrWatson
@quickactivate "JpsiJpsi"

using JLD2
using Plots
theme(:wong)

function get_fit_results(filename)
   f = jldopen(filename)
   JLD2.read(f, "H1_t"), JLD2.read(f, "TS")
end

wl = wload(datadir("simulations", "groups_cross_testing", "Natt=5_Nev=500_Nsampl=500_gen_group=1.jld2"));
fr_0p = wl["fit_results"];

let bins=range(350,500,length=100)
    plot(size=(500,350), xlab="LLH(M{ĥ})", ylab="# sample entries")
    stephist!(map(v->v[1].LLH, fr_0p), bins=bins, lab="fit by group-I")
    stephist!(map(v->v[2].LLH, fr_0p), bins=bins, lab="fit by group-II")
    stephist!(map(v->v[3].LLH, fr_0p), bins=bins, lab="fit by group-III")
    stephist!(map(v->v[4].LLH, fr_0p), bins=bins, lab="fit by group-IV")
end
savefig(plotsdir("llh_testing_higgs.pdf"))

_, fr_1p = get_fit_results(datadir("simulations", "groups_cross_testing", "results_J3_Nev500_Natt5_Nsamples501.jld2"));
let bins=range(370,420,length=100)
    plot(size=(500,350))
    stephist!(map(v->v[1].LLH, fr_1p), bins=bins, lab="group-I")
    stephist!(map(v->v[2].LLH, fr_1p), bins=bins, lab="group-II")
    stephist!(map(v->v[3].LLH, fr_1p), bins=bins, lab="group-III")
    stephist!(map(v->v[4].LLH, fr_1p), bins=bins, lab="group-IV")
end

let bins=range(-30,70,length=50)
    plot(size=(500,350), xlab="TS(0+/1-)", ylab="# sample entries")
    stephist!(-map(v->v[1].LLH - v[3].LLH, fr_0p), bins=bins, lab="0+ generated", seriescolor=1)
    stephist!(-map(v->v[1].LLH - v[3].LLH, fr_1p), bins=bins, lab="1- generated", seriescolor=3)
end
savefig(plotsdir("TS_0p_vs_1m.pdf"))
