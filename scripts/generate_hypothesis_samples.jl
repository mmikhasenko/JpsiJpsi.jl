using JpsiJpsi
using Parameters
using LinearAlgebra
using Plots
using Optim
theme(:wong)

H1_t = randH(1)

@time TS = [let Nev = 500, Natt = 3
    S1_t = sample(Nev; H = H1_t)
    fits = [[fit_sample(S1_t; J = J) for _ in 1:Natt] for J in 0:2]
    best_fits = [f[findmin(getproperty.(f, :LLH))[2]] for f in fits]
    #
    best_fits
end for _ in 1:40*6]

pyplot()
let
    plot(xlab="LLH", ylab="")
    ps = [histogram!(-getproperty.(getindex.(TS,i),:LLH), bins=range(-450,-350,length=100), lab="J=$(i-1)") for i in 1:3]
    plot!(title="fit of J=1 sample", size=(500,350))
end
savefig(joinpath("plots", "llh_of_fit_of_J1.pdf"))

let
    plot(xlab="LLH(S|1)-LLH(S|2)", ylab="LLH(S|2)-LLH(S|3)",
        title="Hyphothesis Testing", size=(500,350))
    hline!([0.0], lab="", l=(:black,0.5))
    vline!([0.0], lab="", l=(:black,0.5))
    scatter!(-(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,2),:LLH)),
            -(getproperty.(getindex.(TS,2),:LLH)-getproperty.(getindex.(TS,3),:LLH)), seriescolor=2,
        lab="J=1")
end
savefig(joinpath("plots", "dLLH_fit_of_J1.pdf"))
#
# histogram(getproperty.(getindex.(TS,2),:LLH)-getproperty.(getindex.(TS,3),:LLH))
# histogram(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,3),:LLH))
# histogram(getproperty.(getindex.(TS,1),:LLH)-getproperty.(getindex.(TS,2),:LLH))
