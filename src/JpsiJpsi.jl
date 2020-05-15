module JpsiJpsi

using PartialWaveFunctions
using Parameters
using Optim

export Hλλ
export nJHs
export LS_for_J0, LS_for_J1, LS_for_J2
include("helicity_couplings.jl")

export xR2vars, vars2xR
export I4μ, I4μ_LS, I4μ_intϕ
export I4K, I4K_LS
include("angular_functions.jl")

export randH
export contract, ×, intensity
export fit_matrix
export algebraic_inversion_matrix
include("fit_of_couplings_matrices.jl")

export randu, randvars
export sample, fit_sample
include("sampling.jl")


end # module
