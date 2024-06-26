"""
    randu() = 2*rand()-1
"""
randu() = 2 * rand() - 1

"""
    randvars() = (cosθ1=randu(), cosθ2=randu(), ϕ=π*randu())
"""
randvars() = (cosθ1=randu(), cosθ2=randu(), ϕ=π * randu())

"""
    sample(Nev; H = error("couplings!"), intensity=I4μ)

produced random sample using hit and miss.
"""
function sample(Nev; H=error("couplings!"), intensity=I4μ)
    f(x) = -intensity(xR2vars(x); H=H)
    search_min(init) = Optim.minimizer(Optim.optimize(f, init, BFGS(),
        Optim.Options(show_trace=false); autodiff=:forwarddiff))
    Imins = f.([search_min(vars2xR(randvars())) for _ in 1:10])
    Imax = -min(Imins...)
    # atan(vars_of_min[1])*2/π, atan(vars_of_min[1])*2/π, atan(vars_of_min[1])*2
    function rand_ev(Imax)
        vars = randvars()
        w = intensity(vars; H=H)
        w > Imax && error("something is wrong: w=$w > Imax=$Imax")
        w < Imax * rand() && return rand_ev(Imax)
        return vars
    end
    return [rand_ev(Imax) for _ in 1:Nev]
end

"""
    fit_sample(S; ng=error("give number of category"), intensity=I4μ)

fit sample `S` using a model of `intensity` with spin hypothesis of category `ng`.
The function returns LLH and the matrix `H` constructed from the fit parameters.
"""
function fit_sample(S; ng=error("give number of category"), intensity=I4μ)
    n = length(ngHs[ng])
    unfold(p) = p[1:n] + 1im .* p[(n+1):end]
    H0(p) = contract(ngHs[ng], unfold(p))
    Hn(p) = H0(p) ./ sqrt(sum(abs2, H0(p)))
    f(pars) = -sum(log, intensity.(S; H=Hn(pars)))
    init_pars = 2 .* rand(2n) .- 1.0
    found_pars = Optim.minimizer(Optim.optimize(f, init_pars, BFGS(),
        Optim.Options(show_trace=false); autodiff=:forwarddiff))
    return (LLH=f(found_pars), H=Hn(found_pars), couplings=unfold(found_pars))
end
