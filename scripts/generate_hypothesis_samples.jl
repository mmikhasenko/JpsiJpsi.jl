using DrWatson
@quickactivate "JpsiJpsi"

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
using Plots
using Optim
using JLD2
theme(:wong)
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
    fulld[:fit_results] = TS
    return fulld
end

params = Dict(
    :Nev => 10,
    :Natt => 1,
    :Nsampl => 2,
    :gen_group => 1,
    :Hgen => randH(1),
    :ng_fit_by => collect(1:7)
)
let
    f = makesim(params)
    # wsave(datadir("simulations", "groups_cross_testing", savename(params, "bson")), f)
    @tagsave(datadir("simulations", "groups_cross_testing", savename(params, "bson")), f)
end

f = makesim(params)
# wsave(datadir("simulations", "groups_cross_testing", savename(params, "bson")), f)
@tagsave(datadir("simulations", "groups_cross_testing", savename(params, "bson")), f)
