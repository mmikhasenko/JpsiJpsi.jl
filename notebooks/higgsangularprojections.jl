### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 297bebd0-1fcc-11ed-3d6d-1f8101b794ad
# ╠═╡ show_logs = false
begin
	using Pkg
	cd(joinpath(@__DIR__, ".."))
	Pkg.activate(".")
	Pkg.instantiate()
	
	using JpsiJpsi
	using Parameters
	using Plots
	using LaTeXStrings
	using Statistics
	using AlgebraPDF
	using FileIO
end

# ╔═╡ a420f22a-2fb0-480e-8ae0-c425cc15e4d8
md"""
# 1D projections of the Higgs sample
"""

# ╔═╡ b89851a9-befd-4c8c-b97d-76295abba197
theme(:wong, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    xlim=(:auto, :auto), ylim=(:auto, :auto))

# ╔═╡ 0c7cbb44-7ad2-48c8-9e08-4907f6e2730b
begin
	@recipe function f(::Type{Val{:errorhist}}, x, y, z)
	    h = Plots._make_hist((y,), plotattributes[:bins], normed=plotattributes[:normalize], weights=plotattributes[:weights])
	    edge = collect(h.edges[1])
	    centers = Plots._bin_centers(edge)
	    dx = Plots.diff(h.edges[1]) / 2
	    x := centers
	    y := h.weights
	    seriestype := :scatter
	    @series begin
	        seriestype := :xerror
	        xerror --> dx
	        yerror --> sqrt.(h.weights)
	        (centers, h.weights)
	    end
	    @series begin
	        seriestype := :yerror
	        yerror --> sqrt.(h.weights)
	        (centers, h.weights)
	    end
	    ()
	end
	@shorthands errorhist
end

# ╔═╡ b2f5d4fa-54a0-4fd7-b84c-3f7a8e8bd909
const Nev = 500

# ╔═╡ cc20ffad-ba31-4e37-9d36-8ea71cb8ce21
md"""
### Generate default sample
"""

# ╔═╡ 1331e9d2-874c-44ee-be5b-1d00ee2238a0
begin
	datasample_filename = datadir("simulations", "higgs_sample", "cosθ12ϕ.jld2")
	if !isfile(datasample_filename)
	    @time S1_t = sample(Nev; H=H_higgs)
	    save(datasample_filename, Dict(:sample => S1_t))
	end
	S1_t = load(datasample_filename)["sample"]
end

# ╔═╡ 3d43d241-1841-4f8c-b7f9-4765e008fc5f
md"""
### Angular modulation PDFs
"""

# ╔═╡ 4b8af723-9c31-4d44-9430-ce3514fd88e5
Iϕ = pdf((@. (ϕ; p) -> 1.0 + p.β * cos(2ϕ)), p=(β=0.3,), lims=(-π, π))

# ╔═╡ e833f755-aa6b-4e3f-aa43-046fc83bf891
Icθi = pdf((@. (c; p) -> 1.0 + p.ζ * (3 * c^2 - 1) / 2), p=(ζ=0.3,), lims=(-1, 1))

# ╔═╡ 695e8cef-af9d-480e-af59-cbfe225b3279
ϕv = getproperty.(S1_t, :ϕ) ;

# ╔═╡ 310fb3d8-6538-4ce1-9687-353a9e49d2c1
md"""
### Fits to the dafault sample
"""

# ╔═╡ 11e1b8b4-59bf-4cfa-963f-df1a92858e81
# ╠═╡ show_logs = false
frβ = fit(Iϕ, ϕv; init_pars=[0.3])[:measurements]

# ╔═╡ cab218ba-d16d-4ca8-ae77-15d459a0cb51
# ╠═╡ show_logs = false
begin
	cosθ1v = getproperty.(S1_t, :cosθ1)
	cosθ2v = getproperty.(S1_t, :cosθ2)
	fit(Icθi, vcat(cosθ1v, cosθ2v); init_pars=[0.0])[:measurements]
end

# ╔═╡ 4f14a6ba-b0bf-4fca-8101-b06f9dd44bd9
md"""
### Statistica uncertainty via sampling
"""

# ╔═╡ 17cfbf43-9569-47e2-8b75-178b2fc60bf5
begin
	betazeta_filename = datadir("simulations", "higgs_sample", "betazeta.jld2")
	if !isfile(betazeta_filename)
	    @time βζs = [
	        let
	            S = sample(Nev; H=H_higgs)
	            ϕv = getproperty.(S, :ϕ)
	            @unpack β = fit(Iϕ, ϕv; init_pars=[0.3])[:parameters]
	            dv = getproperty.(S, :cosθ1)
	            @unpack ζ = fit(Icθi, dv; init_pars=[0.0])[:parameters]
	            (β=β, ζ=ζ)
	        end for _ in 1:500
	    ]
	    #
	    save(betazeta_filename, Dict(:βζs => βζs))
	end
	@unpack βζs = load(betazeta_filename)
end

# ╔═╡ 55cc3d24-3ce7-4cd7-8642-28c6155fa73f
md"""
## Summaty plots
"""

# ╔═╡ f7fdbb8f-27ae-42a0-8ce2-1ea5ca484f6b
# ╠═╡ show_logs = false
let 
	fs = [(f=ϕ -> Iϕ(ϕ, p=(β=1 / 6,)), lab=L"0^+"),
    (f=ϕ -> Iϕ(ϕ, p=(β=-1 / 4,)), lab=L"0^-"),
    (f=ϕ -> Iϕ(ϕ, p=(β=0.0,)), lab=L"1^-"),
	];
	# 
	Nbins = 10
    Nev = length(S1_t)
	# 
    plot(size=(450, 300), xlab=L"\Delta \phi", ylab=L"\#\,\,\mathrm{events}", xlim=(-π, π), leg=:top)
    [plot!(ϕ -> f(ϕ) * (Nev / Nbins * (2π)), -π, π, lab=lab) for (f, lab) in fs]
    errorhist!(ϕv, bins=range(-π, π, length=Nbins + 1), seriescolor=:black, lab=L"\mathrm{data}")
    plot!(ylim=(0, 70), leg=:bottomright)
    #
    plot!(inset=(1, bbox(0.15, 0.51, 0.3, 0.35)))
    plot!(sp=2, xlab=L"\beta", ylab=L"\#\,\,\mathrm{samples}")
    stephist!(sp=2, getproperty.(βζs, :β), frame=:box, lab="", lc=:black,
        bins=range(-1 / 4, 0.32, length=35),
        guidefonthalign=:center)
    vline!(sp=2, [1 / 6], lab="", seriescolor=1, lw=2)
    vline!(sp=2, [-1 / 4], lab="", seriescolor=2, lw=2)
    vline!(sp=2, [0], lab="", seriescolor=3, lw=2,
        xticks=([-1 / 4, 0, 1 / 6], ["-1/4", "0", "1/6"]))
	# 
	savefig(plotsdir("phi_higgs.pdf"))
	plot!()
end

# ╔═╡ 5f549602-2ea9-4eaf-8677-c9dee686f612
# ╠═╡ show_logs = false
let 
	gs = [(f=cosθ -> Icθi(cosθ, p=(ζ=0,)), lab=L"0^{+}"),
	    (f=cosθ -> Icθi(cosθ, p=(ζ=1 / 2,)), lab=L"0^{-}"),
	    (f=cosθ -> Icθi(cosθ, p=(ζ=-1 / 4,)), lab=L"1^{-}"),
	];
	# 
	Nbins = 10
    Nev = length(S1_t)
    plot(size=(450, 300), xlab=L"\cos\theta_1", ylab=L"\#\,\,\mathrm{events}", xlim=(-1, 1), leg=:top)
    [plot!(cosθ -> f(cosθ) * (Nev / Nbins * (2)), -1, 1, lab=lab) for (f, lab) in gs]
    errorhist!(getproperty.(S1_t, :cosθ1), bins=range(-1, 1, length=Nbins + 1), seriescolor=:black, lab=L"\mathrm{data}")
    plot!(ylim=(0, 75), leg=:bottomright)
    #
    plot!(inset=(1, bbox(0.15, 0.51, 0.3, 0.35)))
    plot!(sp=2, xlab="", ann=(0.1, -14, text(L"\zeta", 10)), ylab=L"\#\,\,\mathrm{samples}")
    stephist!(sp=2, getproperty.(βζs, :ζ), frame=:box, lab="", lc=:black,
        bins=range(-1 / 4 - 0.04, 1 / 2 + 0.04, length=35), bottom_margin=0Plots.PlotMeasures.mm)
    vline!(sp=2, [0], lab="", seriescolor=1, lw=2)
    vline!(sp=2, [1 / 2], lab="", seriescolor=2, lw=2)
    vline!(sp=2, [-1 / 4], lab="", seriescolor=3, lw=2,
        xticks=([-1 / 4, 0, 1 / 2], ["-1/4", "0", "1/2"]))
	# 
	savefig(plotsdir("costheta_higgs.pdf"))
	plot!()
end

# ╔═╡ Cell order:
# ╟─a420f22a-2fb0-480e-8ae0-c425cc15e4d8
# ╠═297bebd0-1fcc-11ed-3d6d-1f8101b794ad
# ╠═b89851a9-befd-4c8c-b97d-76295abba197
# ╠═0c7cbb44-7ad2-48c8-9e08-4907f6e2730b
# ╠═b2f5d4fa-54a0-4fd7-b84c-3f7a8e8bd909
# ╟─cc20ffad-ba31-4e37-9d36-8ea71cb8ce21
# ╠═1331e9d2-874c-44ee-be5b-1d00ee2238a0
# ╟─3d43d241-1841-4f8c-b7f9-4765e008fc5f
# ╠═4b8af723-9c31-4d44-9430-ce3514fd88e5
# ╠═e833f755-aa6b-4e3f-aa43-046fc83bf891
# ╠═695e8cef-af9d-480e-af59-cbfe225b3279
# ╟─310fb3d8-6538-4ce1-9687-353a9e49d2c1
# ╠═11e1b8b4-59bf-4cfa-963f-df1a92858e81
# ╠═cab218ba-d16d-4ca8-ae77-15d459a0cb51
# ╟─4f14a6ba-b0bf-4fca-8101-b06f9dd44bd9
# ╠═17cfbf43-9569-47e2-8b75-178b2fc60bf5
# ╟─55cc3d24-3ce7-4cd7-8642-28c6155fa73f
# ╠═f7fdbb8f-27ae-42a0-8ce2-1ea5ca484f6b
# ╠═5f549602-2ea9-4eaf-8677-c9dee686f612
