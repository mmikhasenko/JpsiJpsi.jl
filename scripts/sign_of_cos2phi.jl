using Plots
using DrWatson

theme(:wong)
pyplot()

const Nd = 500;

Mcos2ϕ(S) = sum(cos(2ϕ) for ϕ in S) / length(S)

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

#                                _|            _|          _|
#  _|      _|      _|    _|_|          _|_|_|  _|_|_|    _|_|_|_|    _|_|_|
#  _|      _|      _|  _|_|_|_|  _|  _|    _|  _|    _|    _|      _|_|
#    _|  _|  _|  _|    _|        _|  _|    _|  _|    _|    _|          _|_|
#      _|      _|        _|_|_|  _|    _|_|_|  _|    _|      _|_|  _|_|_|
#                                          _|
#                                      _|_|

function randMcos2ϕ(J, Nev)
    H = randH(J)
    S = [randvars() for _ in 1:Nev]
    Iv = I4μ.(S; H = H)
    return sum(cos(2ϕ)*i for (ϕ,i) in zip(getproperty.(S, :ϕ), Iv)) / sum(Iv)
end

@time M_for_J0_w = [randMcos2ϕ(0, Nd) for _ in 1:1000];
@time M_for_J1_w = [randMcos2ϕ(1, Nd) for _ in 1:1000];
@time M_for_J2_w = [randMcos2ϕ(2, Nd) for _ in 1:1000];

let bins = range(-0.2,0.2,length=50)
    plot(title="Odd and even spins",
        xlab="M2 moment")
    histogram!(M_for_J0_w, bins=bins, lab="J=0", α=0.8)
    histogram!(M_for_J1_w, bins=bins, lab="J=1", α=0.8)
    histogram!(M_for_J2_w, bins=bins, lab="J=2", α=0.8)
end
savefig(joinpath("plots", "moment_M2phi_Nev=$(Nd).pdf"))

const Nd′ = 1500;
@time M_for_J0_w′ = [randMcos2ϕ(0, Nd′) for _ in 1:1000];
@time M_for_J1_w′ = [randMcos2ϕ(1, Nd′) for _ in 1:1000];
@time M_for_J2_w′ = [randMcos2ϕ(2, Nd′) for _ in 1:1000];

let bins = range(-0.2,0.2,length=50)
    plot(title="Odd and even spins",
        xlab="M2 moment")
    histogram!(M_for_J0_w′, bins=bins, lab="J=0", α=0.8)
    histogram!(M_for_J1_w′, bins=bins, lab="J=1", α=0.8)
    histogram!(M_for_J2_w′, bins=bins, lab="J=2", α=0.8)
end
savefig(joinpath("plots", "moment_M2phi_Nev=$(Nd′).pdf"))
