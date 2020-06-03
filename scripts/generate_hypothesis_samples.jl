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
datadir("simulations", "groups_cross_testing")

# path_to_save = joinpath("C:\\","Users","mikha","Documents","JpsiJpsi_material","hypothesis_testing");
#
# Nev = 500; Natt = 5; Nsampl = 500; Jsample = 1
function makesim(d::Dict)
    @unpack ng_gen, ng_fit_by, Nev, Natt, Nsampl = d
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
    :gen_group => 1
    :Hgen => randH(1),
    :ng_fit_by => collect(1:7),
)
let
    f = makesim(params)
    wsave(datadir("simulations", "groups_cross_testing", savename(params, "bson")), f)
    # @tagsave(datadir("simulations", "groups_cross_testing", savename(params, "bson")), f)
end
#
# fs = readdir(datadir("simulations","groups_cross_testing"))[2]
# l = wload(datadir("simulations", "groups_cross_testing", fs))
# l[:fit_results]

# @save joinpath(path_to_save, "results_J$(ngsample)_Nev$(Nev)_Natt$(Natt)_Nsamples$(Nsampl).jld2") f

#
# pyplot()
# let
#     plot(xlab="LLH", ylab="")
#     ps = [stephist!(getproperty.(getindex.(TS,i),:LLH), bins=range(350,450,length=100), lab="by g=$(i)") for i in 1:4]
#     # ps = [stephist!(getproperty.(getindex.(TS,i),:LLH), bins=range(350,450,length=100), lab="J=$(i)", ls=:dash) for i in 5:7]
#     plot!(title="fit of group-I sample", size=(500,350))
# end
# savefig(joinpath("plots", "llh_of_fit_of_J1.pdf"))
# #
# let
#     plot(xlab="LLH(S|1)-LLH(S|2)", ylab="LLH(S|1)-LLH(S|3)",
#         title="Hyphothesis Testing", size=(500,350))
#     hline!([0.0], lab="", l=(:black,0.5))
#     vline!([0.0], lab="", l=(:black,0.5))
#     scatter!((getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,2),:LLH)),
#              (getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,3),:LLH)), seriescolor=1,
#         lab="group 1 generated")
# end
# savefig(joinpath("plots", "dLLH_fit_of_J1.pdf"))
# #
# # histogram(getproperty.(getindex.(TS,2),:LLH)-getproperty.(getindex.(TS,3),:LLH))
# # histogram(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,3),:LLH))
# # histogram(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,2),:LLH))
