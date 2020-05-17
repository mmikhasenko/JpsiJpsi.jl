using JLD2
using

@macroexpand(@load joinpath(path_to_save, "results_J$(Jsample)_Nev$(Nev)_Natt$(Natt)_Nsamples$(Nsampl).jld2") H1_t TS)

function get_fit_results(;
   Jsample = error("J generated"), Nev = error("N events"),
   Natt = 5 , Nsampl = error("size of sample"))
   #
   f = jldopen(joinpath(path_to_save, "results_J$(Jsample)_Nev$(Nev)_Natt$(Natt)_Nsamples$(Nsampl).jld2"))
   JLD2.read(f, "H1_t"), JLD2.read(f, "TS")
end

_, FR0 = get_fit_results(Jsample = 0, Nev = 500, Nsampl = 500);
_, FR1 = get_fit_results(Jsample = 1, Nev = 500, Nsampl = 500);
_, FR2 = get_fit_results(Jsample = 2, Nev = 500, Nsampl = 500);

pyplot()
let
    ps = [begin
        plot(xlab="Log Likelihood of the best fit", ylab="")
        minmax = extrema(vcat([-getproperty.(getindex.(fr,i),:LLH) for i in 1:3]...))
        [histogram!(-getproperty.(getindex.(fr,i),:LLH), bins=range(minmax...,length=100),
            lab="J=$(i-1)", title="generated J=$(j-1) sample", α=0.7) for i in 1:3]
        plot!()
    end for (j,fr) in enumerate([FR0, FR1, FR2])]
    plot(ps..., layout=grid(3,1), size=(700,1000))
end
savefig(joinpath("plots", "llh_of_fit_of_J012.pdf"))
#
let
    plot(xlab="LLH(S|0)-LLH(S|1)", ylab="LLH(S|1)-LLH(S|2)",
        title="Hyphothesis Testing", size=(500,350))
    hline!([0.0], lab="", l=(:black,0.5))
    vline!([0.0], lab="", l=(:black,0.5))
    [
    scatter!(-(getproperty.(getindex.(fr,1),:LLH)-getproperty.(getindex.(fr,2),:LLH)),
            -(getproperty.(getindex.(fr,2),:LLH)-getproperty.(getindex.(fr,3),:LLH)),
        seriescolor=j, lab="J=$(j-1)") for (j,fr) in enumerate([FR0, FR1, FR2])]
    plot!()
end
savefig(joinpath("plots", "dLLH_fit_of_J012.pdf"))
