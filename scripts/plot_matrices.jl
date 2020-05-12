# plotting
let J = 2
    plot(heatmap.(nHs[J+1], c=:balance, clim=(-1,1))...,
        layout=grid(1,length(nHs[J+1])), size=(200*length(nHs[J+1]),190), colorbar=false,
        xaxis=false, yaxis=false)
end
let J = 1
    plot(heatmap.(nHs[J+1], c=:balance, clim=(-1,1))...,
        layout=grid(1,length(nHs[J+1])), size=(200*length(nHs[J+1]),190), colorbar=false,
        xaxis=false, yaxis=false)
end
let J = 0
    plot(heatmap.(nHs[J+1], c=:balance, clim=(-1,1))...,
        layout=grid(1,length(nHs[J+1])), size=(200*length(nHs[J+1]),190), colorbar=false,
        xaxis=false, yaxis=false)
end
