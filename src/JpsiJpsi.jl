module JpsiJpsi

using PartialWaveFunctions
using Parameters
using Optim

export Hλλ
export nJHs
export LS_for_J0, LS_for_J1, LS_for_J2
include("helicity_couplings.jl")

export A4μ, A4μ_intϕ
export A4K
include("angular_functions.jl")

export contract, ×, intensity
export fit_matrix
export algebraic_inversion_matrix
include("fit_of_couplings_matrices.jl")

end # module
