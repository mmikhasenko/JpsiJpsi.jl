module JpsiJpsi

using PartialWaveFunctions
using Parameters
using Optim
using LinearAlgebra
import LinearAlgebra: ×

export H_higgs
include("constants.jl")

export Hλλ
export gHs, ngHs
export groups, specs
include("helicity_couplings.jl")

export xR2vars, vars2xR
export I4μ
export I4K
include("angular_functions.jl")

export randH
export contract, ×, intensity
export fit_matrix
export algebraic_inversion_matrix
include("fit_of_couplings_matrices.jl")

export randu, randvars
export sample, fit_sample
include("sampling.jl")

export datadir, plotsdir
include("io.jl")

end # module
