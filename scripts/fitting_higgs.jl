using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using JpsiJpsi
using Parameters
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
    :Nev => 500,
    :Natt => 5,
    :Nsampl => 500,
    :gen_group => 1,
    :Hgen => H_higgs,
    :ng_fit_by => collect(1:7),
)

f = makesim(params)
save(datadir("simulations", "groups_cross_testing",
        "Natt=$(f[:Natt])_Nev=$(f[:Nev])_Nsampl=$(f[:Nsampl])_gen_group=$(f[:gen_group]).jld2"),
    f)
