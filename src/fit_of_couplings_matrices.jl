contract(H,pars) = sum(h*p for (h,p) in zip(H, pars))
×(H1,H2) = sum(H1 .* conj(H2))

"""
    algebraic_inversion_matrix(H; J=error("give J"))

finds coefficients in the J-basis that approximate H.
Works in the real space.
"""
function algebraic_inversion_matrix(H; J=error("give J"))
    M = [hi × hj for hi in nJHs[J+1], hj in nJHs[J+1]]
    B = [H × hj for hj in nJHs[J+1]]
    return transpose(B \ M)
end

"""
    fit_matrix(H; J = error("give J"))

finds coefficients in the J-basis that approximate H using a fit.
Aslo works in the complex space.
"""
function fit_matrix(H; J = error("give J"))
    n = length(nJHs[J+1])
    unfold(p) = p[1:n] + 1im .* p[(n+1):end]
    f(pars) = intensity(contract(nJHs[J+1], unfold(pars)) - H)
    init_pars = rand(2n)
    found_pars = Optim.minimizer(Optim.optimize(f, init_pars, BFGS(),
                Optim.Options(show_trace = false); autodiff = :forwarddiff))
    return found_pars
end

"""
    randH(J)

generates a random helicity-coupling matrix with complex coefficients
using the basis of spin J.
"""
function randH(J)
    fixed_pars = 2 .* rand(Complex{Float64},length(nJHs[J+1])) .- (1.0+1im)
    fixed_pars ./= sqrt(sum(abs2, fixed_pars))
    H = contract(nJHs[J+1], fixed_pars)
    return H
end
