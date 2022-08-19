using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using JpsiJpsi
using Parameters: @unpack
using LinearAlgebra
using Plots
using Optim
using FileIO
# 
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto))
#
function makesim(d::Dict)
    @unpack Hgen, ng_fit_by, Nev, Natt, Nsampl = d
    #
    fulld = copy(d)
    #
    @time TS = [
        let
            println("\nitteration $(e)/$(Nsampl)\n")
            #
            S1_t = sample(Nev; H=Hgen)
            fits = [[fit_sample(S1_t; ng=ng) for _ in 1:Natt] for ng in ng_fit_by]
            best_fits = [f[findmin(getproperty.(f, :LLH))[2]] for f in fits]
            #
            best_fits
        end for e in 1:Nsampl
    ]
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

f = makesim(params)
save(datadir("simulations",
        "groups_cross_testing",
        "Natt=$(params[:Natt])_Nev=$(params[:Nev])_Nsampl=$(params[:Nsampl])_gen_group=$(params[:gen_group]).jld2"), "results", f)

