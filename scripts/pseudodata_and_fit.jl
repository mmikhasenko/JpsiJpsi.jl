using JpsiJpsi
using Parameters
using LinearAlgebra
using Plots
using Optim
theme(:wong)

H1_t = randH(1)
@time S1_t = sample(500; H = H1_t)

histogram2d(getproperty.(S1_t,:cosθ1),
            getproperty.(S1_t,:cosθ2), bins=10)

function plot_coscos(H; intensity=I4μ_intϕ)
    cosθ1v = range(-1,1,length=110)
    cosθ2v = range(-1,1,length=100)
    calv = [intensity((cosθ1=cosθ1, cosθ2=cosθ2); H=H) for cosθ2 in cosθ2v, cosθ1 in cosθ1v]
    heatmap(cosθ1v, cosθ2v, calv)
end
plot_coscos(H1_t)

p10_t = fit_sample(S1_t; ng = 1)
p11_t = fit_sample(S1_t; ng = 2)
p12_t = fit_sample(S1_t; ng = 3)

-sum(log, I4μ.(S1_t; H = H1_t))

# plot(
#     heatmap(real.(p12_t.H)),
#     heatmap(imag.(p12_t.H)),
#     heatmap(real.(H1_t)),
#     heatmap(imag.(H1_t)))

plot(plot_coscos(H1_t),
     plot_coscos(p10_t.H),
     plot_coscos(p11_t.H),
     plot_coscos(p12_t.H)
)
let l = [p10_t.LLH, p11_t.LLH, p12_t.LLH]
    bar(l .- min(l...))
end
