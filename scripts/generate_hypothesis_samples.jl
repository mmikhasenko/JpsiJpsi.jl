using JpsiJpsi
using Parameters
using LinearAlgebra
using Plots
using Optim
using JLD2
theme(:wong)
#
path_to_save = joinpath("C:\\","Users","mikha","Documents","JpsiJpsi_material","hypothesis_testing");
#
# Nev = 500; Natt = 5; Nsampl = 500; Jsample = 1
const Nev = 500
const Natt = 5
const Nsampl = 501
const ngsample = 1
#
H1_t = randH(ngsample) # J = 1
#
@time TS = [let
    println("\nitteration $(e)/$(Nsampl)\n")
    #
    S1_t = sample(Nev; H = H1_t)
    fits = [[fit_sample(S1_t; ng = ng) for _ in 1:Natt] for ng in 1:7]
    best_fits = [f[findmin(getproperty.(f, :LLH))[2]] for f in fits]
    #
    best_fits
end for e in 1:Nsampl]
#
@save joinpath(path_to_save, "results_J$(ngsample)_Nev$(Nev)_Natt$(Natt)_Nsamples$(Nsampl).jld2") H1_t TS
#
# pyplot()
let
    plot(xlab="LLH", ylab="")
    ps = [stephist!(getproperty.(getindex.(TS,i),:LLH), bins=range(350,450,length=100), lab="by g=$(i)") for i in 1:4]
    # ps = [stephist!(getproperty.(getindex.(TS,i),:LLH), bins=range(350,450,length=100), lab="J=$(i)", ls=:dash) for i in 5:7]
    plot!(title="fit of group-I sample", size=(500,350))
end
savefig(joinpath("plots", "llh_of_fit_of_J1.pdf"))
#
let
    plot(xlab="LLH(S|1)-LLH(S|2)", ylab="LLH(S|1)-LLH(S|3)",
        title="Hyphothesis Testing", size=(500,350))
    hline!([0.0], lab="", l=(:black,0.5))
    vline!([0.0], lab="", l=(:black,0.5))
    scatter!((getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,2),:LLH)),
             (getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,3),:LLH)), seriescolor=1,
        lab="group 1 generated")
end
# savefig(joinpath("plots", "dLLH_fit_of_J1.pdf"))
# #
# # histogram(getproperty.(getindex.(TS,2),:LLH)-getproperty.(getindex.(TS,3),:LLH))
# # histogram(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,3),:LLH))
# # histogram(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,2),:LLH))
