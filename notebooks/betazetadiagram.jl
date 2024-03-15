### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 3f2ba340-1fcb-11ed-0fb1-b714daeb5805
# ╠═╡ show_logs = false
begin
	using Pkg
	cd(joinpath(@__DIR__, ".."))
	Pkg.activate(".")
	Pkg.instantiate()
	
	using JpsiJpsi
	using Parameters
	using LinearAlgebra
	using Plots
	using LaTeXStrings
	using DelimitedFiles
end

# ╔═╡ 205703e3-0694-430d-98ec-38f17e97e059
md"""
# $\zeta/\beta$ correlation diagram
"""

# ╔═╡ 5e84000b-d159-4eea-b43a-8c862ef35333
theme(:wong, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    xlim=(:auto, :auto), ylim=(:auto, :auto))

# ╔═╡ 1afe29d7-226c-4c85-8a31-7256d118bc02
let
    plot(size=(450, 300),
        xlim=(-1.5, 2.5), ylim=(-1.2, 1.2),
        frame=:origin, leg=:bottomright)
    annotate!([(0.1, 1.1, text(L"\beta", 10, :left)), (2.3, -0.05, text(L"\zeta", 10, :top))])
    plot!(ξ -> (2 - ξ) / 3, -1, 2, lab=L"\mathrm{group}\,\,I", fill_between=0, seriescolor=1)
    plot!(ξ -> (2 - ξ) / 3, -1, 2, lab="", lw=2, seriescolor=1)
    hline!([0.0], lab="", lc=:black)
    vline!([0.0], lab="", lc=:black)
    plot!(ξ -> -(1 - 2ξ) / 3, -1, 1 / 2, lab=L"\mathrm{group}\,\,I\!I", lw=6, seriescolor=2)
    scatter!([1 / 2], [0], lab=L"\mathrm{group}\,\,I\!I\!I", ms=13, markerstrokewidth=0, markerstrokecolor=3, lc=:white, seriescolor=3)
    plot!(ξ -> 0, -1, 1 / 2, lab=L"\mathrm{group}\,\,I\!V", lw=6, seriescolor=4)
    scatter!([1 / 2], [0], lab="", ms=13, markerstrokewidth=0, markerstrokecolor=3, seriescolor=3)
	# 
	savefig(plotsdir("diagram_K.pdf"))
	plot!()
end

# ╔═╡ d89ab7c6-8546-470d-a2bb-68b150fe960a
let
    plot(size=(450, 300),
        xlim=(-1.3, 0.8), ylim=(-0.3, 0.3),
        yticks=([-1 / 4, 0, 1 / 4]),
        frame=:origin, leg=:bottomleft)
    annotate!([(0.05, 0.27, text(L"\beta", 10, :left)), (0.7, -0.01, text(L"\zeta", 10, :top))])
    plot!(ξ -> (1 + ξ) / 6, -1, 1 / 2, lab=L"\mathrm{group}\,\,I", fill_between=0, seriescolor=1)
    plot!(ξ -> (1 + ξ) / 6, -1, 1 / 2, lab="", lw=2, seriescolor=1)
    hline!([0.0], lab="", lc=:black)
    vline!([0.0], lab="", lc=:black)
    # plot!(ξ->-(1/2+2ξ)/6, -1/4, 1/2, lab="", fill_between=0, α=0.5, seriescolor=2)
    plot!(ξ -> -(1 / 2 + 2ξ) / 6, -1 / 4, 1 / 2, lab=L"\mathrm{group}\,\,I\!I", lw=6, seriescolor=2)
    scatter!([-1 / 4], [0], lab=L"\mathrm{group}\,\,I\!I\!I", ms=13, markerstrokewidth=0, markerstrokecolor=3, lw=0, seriescolor=3)
    plot!(ξ -> 0, -1 / 4, 1 / 2, lab=L"\mathrm{group}\,\,I\!V", lw=6, seriescolor=4)
    scatter!([-1 / 4], [0], lab="", ms=13, markerstrokewidth=0, markerstrokecolor=3, lw=0, seriescolor=3)
	# 
	savefig(plotsdir("diagram_mu.pdf"))
	plot!()
end

# ╔═╡ Cell order:
# ╟─205703e3-0694-430d-98ec-38f17e97e059
# ╠═3f2ba340-1fcb-11ed-0fb1-b714daeb5805
# ╠═5e84000b-d159-4eea-b43a-8c862ef35333
# ╠═1afe29d7-226c-4c85-8a31-7256d118bc02
# ╠═d89ab7c6-8546-470d-a2bb-68b150fe960a
