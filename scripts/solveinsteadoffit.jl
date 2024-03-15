using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()

using JpsiJpsi
using Parameters
using LinearAlgebra
using Plots
using Optim

theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto, :auto), ylim=(:auto, :auto))

#
function solve_mismatch_real(; ng=1, ng′=2)
    H = randH(ng)
    found_pars = algebraic_inversion_matrix(H; ng=ng′)
    mismatch = intensity(contract(ngHs[ng′], transpose(found_pars)) - H)
    return mismatch
end

function fit_mismatch(; ng=1, ng′=2)
    H = randH(ng)
    found_pars = fit_matrix(H; ng=ng′)
    n = length(ngHs[ng′])
    unfold(p) = p[1:n] + 1im .* p[(n+1):end]
    mismatch = intensity(contract(ngHs[ng′], transpose(unfold(found_pars))) - H)
    return mismatch
end
#
histogram([fit_mismatch(; ng=1, ng′=5) for _ in 1:1000], bins=range(0, 2, length=100))
histogram([fit_mismatch(; ng=2, ng′=6) for _ in 1:1000], bins=range(0, 2, length=100))
histogram([fit_mismatch(; ng=4, ng′=7) for _ in 1:1000], bins=range(0, 2, length=100))
# 
