# JpsiJpsi

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12583722.svg)](https://doi.org/10.5281/zenodo.12583722)
[![inspire reference](https://img.shields.io/badge/article-JHEP_0424-green)](https://inspirehep.net/literature/1806437)


The repository collects material for exploratory studies of J/ψ (μ⁺μ⁻) J/ψ (μ⁺μ⁻) system.

## Repositories

- [JpsiJpsi](https://github.com/mmikhasenko/JpsiJpsi.jl)

## Repository content

Paper plots in [`scripts`](scripts/) folder:

- `fitdata.jl`: fit CMS data with incoherent model

ZZ analysis:

- `higgsgenerateandfit.jl`: generate and fit samples
- `higgshypothesestest.jl`: paper figures
- `higgsangularprojections.jl`: theta, phi distributions for Higgs

Additional exercises:

- `theta1theta2correlations.jl`: cosθ₁-cosθ₂ correlation
- `cos2phimoments.jl`: check beta-sign using moments for randH

### Pluto notebooks:

- `betazetadiagram.jl`: beta-zeta diagrams
- `massdepbetazeta.jl`: mass-beta-zeta plots for different mass hypotheses

To run the scripts just open pluto browser.

```julia
julia> using Pluto; Pluto.run()
```
