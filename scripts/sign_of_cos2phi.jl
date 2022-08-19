using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using LinearAlgebra

using Plots
theme(:wong)
# pyplot()

const Nd = 500;

Mcos2ϕ(S) = sum(cos.(2 .* S)) / length(S)

S1_t = [Mcos2ϕ(getproperty.(sample(Nd; H=Diagonal(fill(1 / sqrt(3), 3))), :ϕ)) for _ in 1:100]
S2_t = [Mcos2ϕ(getproperty.(sample(Nd; H=Diagonal([1 / √2, 0, -1 / √2])), :ϕ)) for _ in 1:100]
S3_t = [Mcos2ϕ(getproperty.(sample(Nd; H=ngHs[3][1]), :ϕ)) for _ in 1:100]

let bins = range(-0.2, 0.2, length=50)
    plot(title="Odd and even spins",
        xlab="M2 moment")
    histogram!(S1_t, bins=bins, lab="group-I", α=0.8)
    histogram!(S2_t, bins=bins, lab="group-II", α=0.8)
    histogram!(S3_t, bins=bins, lab="group-III", α=0.8)
end
savefig(joinpath("plots", "moment_M2phi_Nev=500.pdf"))

#                                _|            _|          _|
#  _|      _|      _|    _|_|          _|_|_|  _|_|_|    _|_|_|_|    _|_|_|
#  _|      _|      _|  _|_|_|_|  _|  _|    _|  _|    _|    _|      _|_|
#    _|  _|  _|  _|    _|        _|  _|    _|  _|    _|    _|          _|_|
#      _|      _|        _|_|_|  _|    _|_|_|  _|    _|      _|_|  _|_|_|
#                                          _|
#                                      _|_|

function randMcos2ϕ(ng, Nev)
    H = randH(ng)
    S = [randvars() for _ in 1:Nev]
    Iv = I4μ.(S; H=H)
    return sum(cos(2ϕ) * i for (ϕ, i) in zip(getproperty.(S, :ϕ), Iv)) / sum(Iv)
end

@time M_for_ng1_w = [randMcos2ϕ(1, Nd) for _ in 1:1000];
@time M_for_ng2_w = [randMcos2ϕ(2, Nd) for _ in 1:1000];
@time M_for_ng3_w = [randMcos2ϕ(3, Nd) for _ in 1:1000];

let bins = range(-0.2, 0.2, length=50)
    plot(title="Odd and even spins",
        xlab="M2 moment")
    histogram!(M_for_ng1_w, bins=bins, lab="group-I", α=0.8)
    histogram!(M_for_ng2_w, bins=bins, lab="group-II", α=0.8)
    histogram!(M_for_ng3_w, bins=bins, lab="group-III", α=0.8)
end
savefig(joinpath("plots", "moment_M2phi_Nev=$(Nd).pdf"))

const Nd′ = 1500;
@time M_for_ng1_w′ = [randMcos2ϕ(0, Nd′) for _ in 1:1000];
@time M_for_ng2_w′ = [randMcos2ϕ(1, Nd′) for _ in 1:1000];
@time M_for_ng3_w′ = [randMcos2ϕ(2, Nd′) for _ in 1:1000];

let bins = range(-0.2, 0.2, length=50)
    plot(title="Odd and even spins",
        xlab="M2 moment")
    histogram!(M_for_ng1_w′, bins=bins, lab="group-I", α=0.8)
    histogram!(M_for_ng2_w′, bins=bins, lab="group-II", α=0.8)
    histogram!(M_for_ng3_w′, bins=bins, lab="group-III", α=0.8)
end
savefig(joinpath("plots", "moment_M2phi_Nev=$(Nd′).pdf"))
