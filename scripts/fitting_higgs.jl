using DrWatson
@quickactivate "JpsiJpsi"

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
#
function makesim(d::Dict)
    @unpack Hgen, ng_fit_by, Nev, Natt, Nsampl = d
    #
    fulld = copy(d)
    #
    @time TS = [let
        println("\nitteration $(e)/$(Nsampl)\n")
        #
        S1_t = sample(Nev; H = Hgen)
        fits = [[fit_sample(S1_t; ng = ng) for _ in 1:Natt] for ng in ng_fit_by]
        best_fits = [f[findmin(getproperty.(f, :LLH))[2]] for f in fits]
        #
        best_fits
    end for e in 1:Nsampl]
    #
    fulld["fit_results"] = TS
    return fulld
end

params = Dict(
    "Nev" => 500,
    "Natt" => 5,
    "Nsampl" => 500,
    "gen_group" => 1,
    "Hgen" => Diagonal(fill(2.0,3)),
    "ng_fit_by" => collect(1:7),
)

f = makesim(params)
wsave(datadir("simulations", "groups_cross_testing", savename(params, "jld2")), f)

let bins=range(350,500,length=100)
    plot(size=(500,350))
    stephist!(map(v->v[1].LLH, f["fit_results"]), bins=bins, lab="group-I")
    stephist!(map(v->v[2].LLH, f["fit_results"]), bins=bins, lab="group-II")
    stephist!(map(v->v[3].LLH, f["fit_results"]), bins=bins, lab="group-III")
    stephist!(map(v->v[4].LLH, f["fit_results"]), bins=bins, lab="group-IV")
end
savefig(plotsdir("llh_testing_higgs.pdf"))
