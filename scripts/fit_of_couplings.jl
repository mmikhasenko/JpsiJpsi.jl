using JpsiJpsi
using Parameters
using LinearAlgebra
using Plots
using Optim
theme(:wong)

import JpsiJpsi: ×

# orthogonality
[hi × hj for hi in ngHs[1], hj in ngHs[3]]
[hi × hj for hi in ngHs[1], hj in ngHs[2]]
[hi × hj for hi in ngHs[3], hj in ngHs[2]]
[sum(hi × hj for hi in ngHs[i], hj in ngHs[j]) for i = 1:7, j = 1:7]

#
function solve_mismatch_real(; ng=1, ng′=2)
    H = randH(ng)
    found_pars = algebraic_inversion_matrix(H; ng=ng′)
    mismatch = intensity(contract(ngHs[ng′], transpose(found_pars)) - H)
end

# histogram(log10.([solve_mismatch_real(; ng = 1, ng′ = 3) for _ in 1:10000]), bins=100, xlab="log10(Delta)")
#
function fit_mismatch(; ng=1, ng′=2)
    H = randH(ng)
    found_pars = fit_matrix(H; ng=ng′)
    n = length(ngHs[ng′])
    unfold(p) = p[1:n] + 1im .* p[(n+1):end]
    mismatch = intensity(contract(ngHs[ng′], transpose(unfold(found_pars))) - H)
end
#
histogram([fit_mismatch(; ng=1, ng′=5) for _ in 1:1000], bins=range(0, 2, length=100))
histogram([fit_mismatch(; ng=2, ng′=6) for _ in 1:1000], bins=range(0, 2, length=100))
histogram([fit_mismatch(; ng=4, ng′=7) for _ in 1:1000], bins=range(0, 2, length=100))
# 
