using JpsiJpsi

using Plots
using LaTeXStrings
pyplot()

let (L, S) = (0, 0), J = 0
    cv = range(-1, 1, length=100)
    v = [A4μ_intϕ((cosθ1=c1, cosθ2=c2); LS=(L, S), J=J) for (c1, c2) in Iterators.product(cv, cv)]
    v ./= sum(v) / (100 / 2)^2
    plot(size=(450, 400),
        xlab=L"\cos\,\theta_1", ylab=L"\cos\,\theta_2", title="J = $J: (L=$L, S=$S)")
    heatmap!(cv, cv, v, c=:viridis)
end

let (L, S) = (0, 0), J = 0
    cv = range(-1, 1, length=100)
    v = [A4μ_intϕ((cosθ1=c1, cosθ2=c2); LS=(L, S), J=J) for (c1, c2) in Iterators.product(cv, cv)]
    v ./= sum(v) / (100 / 2)^2
    plot(size=(450, 400),
        xlab=L"\cos\,\theta_1", ylab=L"\cos\,\theta_2", title="J = $J: (L=$L, S=$S)")
    heatmap!(cv, cv, v, c=:viridis)
end

for (J, LSs) in (0 => LS_for_J0, 1 => LS_for_J1, 2 => LS_for_J2)
    for (L, S) in LSs
        cv = range(-1, 1, length=100)
        v = [A4μ_intϕ((cosθ1=c1, cosθ2=c2); LS=(L, S), J=J) for (c1, c2) in Iterators.product(cv, cv)]
        v ./= sum(v) / (100 / 2)^2
        plot(size=(450, 400),
            xlab=L"\cos\,\theta_1", ylab=L"\cos\,\theta_2", title="J = $J: (L=$L, S=$S)")
        heatmap!(cv, cv, v, c=:viridis, colorbar=false)
        savefig(joinpath("plots", "map_JLS_$(J)$(L)$(S).pdf"))
    end
end
