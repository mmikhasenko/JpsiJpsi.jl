using Plots
theme(:wong)
pyplot()

H1_t = randH(1) # J = 1
S1_t = sample(500; H = H1_t);
histogram(getproperty.(S1_t, :ϕ), bins=50)

const Nd = 500;

Mcos2ϕ(S) = sum(cos(2ϕ) for ϕ in S)

@time M_for_J0 = [Mcos2ϕ(getproperty.(sample(Nd; H = randH(0)), :ϕ)) for _ in 1:1000];
@time M_for_J1 = [Mcos2ϕ(getproperty.(sample(Nd; H = randH(1)), :ϕ)) for _ in 1:1000];
@time M_for_J2 = [Mcos2ϕ(getproperty.(sample(Nd; H = randH(2)), :ϕ)) for _ in 1:1000];

histogram(M_for_J0, bins=100)
histogram(M_for_J1, bins=100)
histogram(M_for_J2, bins=100)

let bins = range(-100,100,length=50)
    plot(title="Odd and even spins",
        xlab="M2 moment")
    histogram!(M_for_J0, bins=bins, lab="J=0", α=0.8)
    histogram!(M_for_J1, bins=bins, lab="J=1", α=0.8)
    histogram!(M_for_J2, bins=bins, lab="J=2", α=0.8)
end
savefig(joinpath("plots", "moment_M2phi_Nev=500.pdf"))
