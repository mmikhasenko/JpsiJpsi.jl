# JpsiJpsi

The repository collects material for exploratory studies of J/ψ (μ⁺μ⁻) J/ψ (μ⁺μ⁻) system.

## Repositories

- [JpsiJpsi]()
- [cep-phiphi-angular-analysis](https://gitlab.cern.ch/LHCb-QEE/cep-phiphi-angular-analysis)

## Mattermost

- [X2VV](https://mattermost.web.cern.ch/lhcb/channels/x2vv)

## Overleaf

- [link to edit](https://www.overleaf.com/2812294971bkkkbqgtcrtg#062367)

## Preview

Plots and computations in [GitLab Pages](http://h2vv.docs.cern.ch)!

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
