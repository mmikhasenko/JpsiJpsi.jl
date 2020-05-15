using JpsiJpsi
using Parameters
using LinearAlgebra
using Plots
using Optim
theme(:wong)

import JpsiJpsi: ×

# orthogonality
[hi × hj for hi in nJHs[1], hj in nJHs[3]]
[hi × hj for hi in nJHs[1], hj in nJHs[2]]
[hi × hj for hi in nJHs[3], hj in nJHs[2]]
#
function solve_mismatch_real(; J = 0, J′ = 2)
    H = randH(J)
    found_pars = algebraic_inversion_matrix(H; J=J′)
    mismatch = intensity(contract(nJHs[J′+1], transpose(found_pars)) - H)
end

histogram(log10.([solve_mismatch_real(; J=2, J′=0) for _ in 1:10000]), bins=100, xlab="log10(Delta)")
#
function fit_mismatch(; J = 2, J′ = 0)
    H = randH(J)
    found_pars = fit_matrix(H; J=J′)
    n = length(nJHs[J′+1])
    unfold(p) = p[1:n] + 1im .* p[(n+1):end]
    mismatch = intensity(contract(nJHs[J′+1], transpose(unfold(found_pars))) - H)
end
#
histogram([fit_mismatch(; J=2, J′=0) for _ in 1:10000], bins=100)
