using JLD2

path_to_save = joinpath("C:\\", "Users", "mikha", "Documents", "JpsiJpsi_material", "hypothesis_testing");

function get_fit_results(;
    Jsample=error("J generated"), Nev=error("N events"),
    Natt=5, Nsampl=error("size of sample"))
    #
    f = jldopen(joinpath(path_to_save, "results_J$(Jsample)_Nev$(Nev)_Natt$(Natt)_Nsamples$(Nsampl).jld2"))
    JLD2.read(f, "H1_t"), JLD2.read(f, "TS")
end

_, FR3 = get_fit_results(Jsample=3, Nev=500, Nsampl=501);
# _, FR1 = get_fit_results(Jsample = 1, Nev = 500, Nsampl = 500);
# _, FR2 = get_fit_results(Jsample = 2, Nev = 500, Nsampl = 500);

FR3

pyplot()
let
    ps = [
        begin
            plot(xlab="Log Likelihood of the best fit", ylab="")
            minmax = extrema(vcat([getproperty.(getindex.(fr, i), :LLH) for i in 1:4]...))
            [histogram!(getproperty.(getindex.(fr, i), :LLH), bins=range(minmax..., length=100),
                lab="group-" * ["I", "II", "III", "IV"][i], title="generated group-III sample", α=0.7) for i in 1:4]
            plot!()
        end for (j, fr) in enumerate([FR3])
    ]
    plot(ps..., size=(500, 350))
end
savefig(joinpath("plots", "llh_of_fit_of_J012.pdf"))

let iΔ_x = (3, 2), iΔ_y = (3, 1)
    plot(xlab="LLH(S|$(iΔ_x[1]))-LLH(S|$(iΔ_x[2]))", ylab="LLH(S|$(iΔ_y[1]))-LLH(S|$(iΔ_y[2 ]))",
        title="Hyphothesis Testing", size=(500, 350))
    hline!([0.0], lab="", l=(:black, 0.5))
    vline!([0.0], lab="", l=(:black, 0.5))
    scatter!((getproperty.(getindex.(FR3, iΔ_x[1]), :LLH) - getproperty.(getindex.(FR3, iΔ_x[2]), :LLH)),
        (getproperty.(getindex.(FR3, iΔ_y[1]), :LLH) - getproperty.(getindex.(FR3, iΔ_y[2]), :LLH)), seriescolor=1,
        lab="group 1 generated")
end


#
let
    plot(xlab="LLH(S|0)-LLH(S|1)", ylab="LLH(S|1)-LLH(S|2)",
        title="Hyphothesis Testing", size=(500, 350))
    hline!([0.0], lab="", l=(:black, 0.5))
    vline!([0.0], lab="", l=(:black, 0.5))
    [
        scatter!(-(getproperty.(getindex.(fr, 1), :LLH) - getproperty.(getindex.(fr, 2), :LLH)),
            -(getproperty.(getindex.(fr, 2), :LLH) - getproperty.(getindex.(fr, 3), :LLH)),
            seriescolor=j, lab="J=$(j-1)") for (j, fr) in enumerate([FR0, FR1, FR2])]
    plot!()
end
savefig(joinpath("plots", "dLLH_fit_of_J012.pdf"))
