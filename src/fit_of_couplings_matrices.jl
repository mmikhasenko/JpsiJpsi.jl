contract(H,pars) = sum(h*p for (h,p) in zip(H, pars))
×(H1,H2) = sum(H1 .* conj(H2))

"""
    algebraic_inversion_matrix(H; ng=error("give number of category"))

finds coefficients in the group basis indexed with `ng` that approximate H.
Works in the real space.
"""
function algebraic_inversion_matrix(H; ng=error("give number of category"))
    M = [hi × hj for hi in ngHs[ng], hj in ngHs[ng]]
    B = [H × hj for hj in ngHs[ng]]
    return transpose(B \ M)
end

"""
    fit_matrix(H; ng=error("give number of category"))

finds coefficients in the group basis indexed with `ng` that approximate H using a fit.
Aslo works in the complex space.
"""
function fit_matrix(H; ng=error("give number of category"))
    n = length(ngHs[ng])
    unfold(p) = p[1:n] + 1im .* p[(n+1):end]
    f(pars) = intensity(contract(ngHs[ng], unfold(pars)) - H)
    init_pars = rand(2n)
    found_pars = Optim.minimizer(Optim.optimize(f, init_pars, BFGS(),
                Optim.Options(show_trace = false); autodiff = :forwarddiff))
    return found_pars
end

"""
    randH(ng)

generates a random helicity-coupling matrix with complex coefficients
using the basis of groups indexed by ng.
"""
function randH(ng)
    fixed_pars = 2 .* rand(Complex{Float64},length(ngHs[ng])) .- (1.0+1im)
    fixed_pars ./= sqrt(sum(abs2, fixed_pars))
    H = contract(ngHs[ng], fixed_pars)
    return H
end
