### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 05ff1b20-0756-11ed-12b2-3118481be1cd
begin
	using SymPy
	using LinearAlgebra
	# 
	using Parameters
	
	import PyCall
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	import_from(sympy.physics.quantum.spin, (:WignerD, :wignerd), typ=:Any)
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	import_from(sympy.physics.wigner)
end

# ╔═╡ c55fcb81-5248-4557-909b-cc306efe89a5
Wignerd(J,λ1,λ2,θ) = WignerD(Sym(J),λ1,λ2,0,θ,0).doit()

# ╔═╡ cdf9201c-e88f-43ee-9bb4-f847f73439ea
md"""
## Decay Amplitude
"""

# ╔═╡ d989eab6-53e8-49d5-9c2b-ecb59269f308
function I4P(vars, ξrange; H)
    @unpack θ1,ϕ,θ2 = vars
	nξ = length(ξrange)
    return simplify(9/Sym(nξ^2)*sum(
        (λ1-λ2 == λ1′-λ2′ ? 1 : 0) *
            Wignerd(1,λ1 ,ξ1,θ1) * Wignerd(1,λ2 ,ξ2,θ2) * (isodd(1-λ2) ? -1 : 1) *
            Wignerd(1,λ1′,ξ1,θ1) * Wignerd(1,λ2′,ξ2,θ2) * (isodd(1-λ2′) ? -1 : 1) *
            (cos((λ1-λ1′)*ϕ)+1im*sin((λ1-λ1′)*ϕ)) *
                 H[λ1+2,λ2+2] *
            conj(H[λ1′+2,λ2′+2])
            for ξ1 in ξrange, ξ2 in ξrange,
                λ1  in -1:1, λ2  in -1:1,
                λ1′ in -1:1, λ2′ in -1:1))
end

# ╔═╡ ae1fb352-5826-480f-a503-e1b6951e0126
a,b,c,d,ϵ,P,s = @vars a b c d ϵ P s real=true

# ╔═╡ 71fabb02-3187-4de9-982b-77c74ea01727
Hg = [b   a     c;
      s*a   d ϵ*s*a;
      s*c ϵ*a   ϵ*b];

# ╔═╡ d50197fe-28b4-440b-8d26-5ae98995e5fe
normalization_relation = sum(Hg.^2).subs(Dict(ϵ^2 => 1, s^2 => 1))

# ╔═╡ c9f00440-2114-4057-854e-8e563f62d52e
md"""
### Two example: leptons and scalars
"""

# ╔═╡ e57abfb3-8a3e-4e51-9039-393110f16af5
begin
	I4μ(vars; H) = I4P(vars, (-1,1); H)
	I4K(vars; H) = I4P(vars, (0,); H)
end

# ╔═╡ 0bac613d-afb8-4cb6-8dbc-226c3ebaff6d
θ1, θ2, ϕ = @vars θ_1 θ_2 ϕ positive=true

# ╔═╡ abc17790-fef5-41d6-8a89-1e53a525d45e
I4μ_full = I4μ((; θ1,ϕ,θ2); H=Hg) ;

# ╔═╡ e361dd95-e168-42fa-843f-815d3ce21f51
I4K_full = I4K((; θ1,ϕ,θ2); H=Hg) ;

# ╔═╡ 0fc4b7db-0589-44de-98aa-18fd76b2126e
md"""
## Expansion in a basis
"""

# ╔═╡ 099913fe-ac20-4d1c-8755-a7ac4d1da78c
const angular_basis = [
    sin(θ1)^2*sin(θ2)^2*sin(ϕ)^2,
    sin(θ1)sin(θ2)cos(θ1)cos(θ2)*cos(ϕ),
    sin(θ1)^2*sin(θ2)^2,
    sin(θ1)^2,
    sin(θ2)^2,
	1
]

# ╔═╡ e7deeab1-a4c1-46d8-8ac1-2da99e1d18d4
begin
	integrateθ(f, θ) = integrate(f*sin(θ),	(θ, 0, sympy.pi)) / 2
	integrateϕ(f) = integrate(f, (ϕ, -sympy.pi, sympy.pi)) / (2sympy.pi)
	integrateθs(f) = f |> e->integrateθ(e, θ1) |> e->integrateθ(e, θ2)
	integrateϕθs(f) = f |> integrateϕ |> integrateθs
end

# ╔═╡ dd13511c-54bd-4c66-99b0-2342da4768d1
norms = integrateϕθs.(angular_basis)

# ╔═╡ 782084b4-6e62-4cad-a81c-cd9d11ab6a1c
function expand_on_basis(ampl, basis)
    lexp = ampl.subs(Dict(
		# 
		cos(θ1)^2 => 1-sin(θ1)^2,
		cos(θ2)^2 => 1-sin(θ2)^2,
		cos(2ϕ) => 1-2*sin(ϕ)^2,
		# 
		ϵ^2 => 1,
		s^2 => 1)) |> sympy.expand_trig |> sympy.expand
	# 
    coeff = Vector{Sym}(undef,length(basis))
    for (i,b) in enumerate(basis[1:end-1]) # exclude basis[end] = 1
		polynomialforb = collect(lexp, b)
		c = polynomialforb.coeff(b, 1) 
        coeff[i] = c
		lexp = lexp - c*b  |> sympy.expand
    end
    coeff[end] = lexp
    return together.(coeff)
end

# ╔═╡ 63442521-9fe1-4319-9f13-82daf14b832f
md"""
### $J/\psi J/\psi \to 4 \mu$
"""

# ╔═╡ a6502774-b01e-49ae-9e1b-5f58373b8525
coeff_μ = expand_on_basis(I4μ_full, angular_basis) ;

# ╔═╡ 3b46e566-c7f2-4892-91b1-59bfd91c91b1
@assert simplify(sum(coeff_μ .* norms)) ==
	sum(Hg.^2).subs(Dict(ϵ^2 => 1, s^2 => 1))

# ╔═╡ d6189855-bed3-4758-b14d-b93849cb04b9
print_coeff_μ = coeff_μ .*  map(x->x==0 ? 1 : x, norms)

# ╔═╡ 710d8bdc-32e8-484a-a698-8ff30c1ad27a
md"""
### $\phi \phi \to 4 K$
"""

# ╔═╡ 02529c4d-79d6-4210-8add-ec2375ac2715
coeff_K = expand_on_basis(I4K_full, angular_basis) ;

# ╔═╡ 05fb0941-85ad-4768-9adb-119f15542c1f
simplify(sum(coeff_K .* norms))

# ╔═╡ 38573cec-0eb5-47eb-b15f-0dc759f6872f
@assert simplify(sum(coeff_K .* norms)) == normalization_relation

# ╔═╡ 3500d280-3fe2-4054-aed5-dd98491216d5
print_coeff_K = together.(coeff_K .*  map(x->x==0 ? 1 : x, norms))

# ╔═╡ 49f63b21-88e5-49b6-b5fc-2708553a9637
md"""
## Asymmetries
"""

# ╔═╡ 87ec64a5-ff9c-4219-87e3-dfeddb628684
md"""
### Beta asymmetry
"""

# ╔═╡ ca2773ba-c781-4786-be3d-d9c0d549a1c8
asymmetryβ_μ = I4μ_full |> integrateθs ;

# ╔═╡ 152ac4eb-7eb5-41f8-b9d7-bffa5e61f4c2
asymmetryβ_K = I4K_full |> integrateθs ;

# ╔═╡ 2016bd95-739a-4227-a4a0-bab46676296c
function getβ(f)
	expr = f |> 
		subs(Dict(ϵ^2=>1,s^2=>1)) |>
		subs(sin(ϕ)^2=>(1-cos(2ϕ))/2) |> simplify
	P = sympy.Poly(expr, cos(2ϕ))
	@assert P.degree() == 1
	N, D = P.coeffs()
	return N/D
end

# ╔═╡ 92d7b42d-3ff9-4516-b695-11c873d429ce
Dict(:βμ => getβ(asymmetryβ_μ))

# ╔═╡ 4650b3c9-c490-4520-a45d-4666b5f31bc8
Dict(:βK => getβ(asymmetryβ_K))

# ╔═╡ df143390-dcd2-4709-b029-a56fc1e0f73c
md"""
### Zeta asymmetry
"""

# ╔═╡ 02aa1aa3-ef47-47be-bbbd-bae4feb73b4b
asymmetryζ_μ = integrateθ(I4μ_full,θ1) |> integrateϕ |> expand;

# ╔═╡ 4b1009c0-5f5a-4373-aaed-0b1a1c9a7bc9
function getζ(f)
	x = Sym("x")
	expr = f |> 
		subs(Dict(ϵ^2=>1,s^2=>1)) |>
		subs(cos(θ2)^2=>1-sin(θ2)^2) |>
		subs(sin(θ2)^2=>1-(2x+1)/3) |> simplify
	P = sympy.Poly(expr, x)
	@assert P.degree() == 1
	N, D = P.coeffs()
	return N/D
end

# ╔═╡ 0770f489-ab7b-4fa3-aeb5-295cff0371b4
Dict(:ζμ => getζ(asymmetryζ_μ))

# ╔═╡ 09f3d176-9529-4fdc-bc74-5578b23dc06f
asymmetryζ_K = integrateθ(I4K_full,θ1) |> integrateϕ ;

# ╔═╡ 250df9a1-a935-4325-b3f1-448a172fec1a
Dict(:ζK => getζ(asymmetryζ_K))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PyCall = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"

[compat]
Parameters = "~0.12.3"
PyCall = "~1.93.1"
SymPy = "~1.1.6"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "ff38036fb7edc903de4e79f32067d8497508616b"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.2"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CommonEq]]
git-tree-sha1 = "d1beba82ceee6dc0fce8cb6b80bf600bbde66381"
uuid = "3709ef60-1bee-4518-9f2f-acd86f176c50"
version = "0.2.0"

[[deps.CommonSolve]]
git-tree-sha1 = "332a332c97c7071600984b3c31d9067e1a4e6e25"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.1"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "924cdca592bc16f14d2f7006754a621735280b74"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.1.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "46a39b9c58749eefb5f2dc1178cb8fab5332b1ab"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.15"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "1fc929f47d7c151c839c5fc1375929766fb8edcc"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.93.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SymPy]]
deps = ["CommonEq", "CommonSolve", "Latexify", "LinearAlgebra", "Markdown", "PyCall", "RecipesBase", "SpecialFunctions"]
git-tree-sha1 = "e1865ba3c44551087a04295ddc40c10edf1b24a0"
uuid = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
version = "1.1.6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
"""

# ╔═╡ Cell order:
# ╠═05ff1b20-0756-11ed-12b2-3118481be1cd
# ╠═c55fcb81-5248-4557-909b-cc306efe89a5
# ╟─cdf9201c-e88f-43ee-9bb4-f847f73439ea
# ╠═d989eab6-53e8-49d5-9c2b-ecb59269f308
# ╠═ae1fb352-5826-480f-a503-e1b6951e0126
# ╠═71fabb02-3187-4de9-982b-77c74ea01727
# ╠═d50197fe-28b4-440b-8d26-5ae98995e5fe
# ╟─c9f00440-2114-4057-854e-8e563f62d52e
# ╠═e57abfb3-8a3e-4e51-9039-393110f16af5
# ╠═0bac613d-afb8-4cb6-8dbc-226c3ebaff6d
# ╠═abc17790-fef5-41d6-8a89-1e53a525d45e
# ╠═e361dd95-e168-42fa-843f-815d3ce21f51
# ╟─0fc4b7db-0589-44de-98aa-18fd76b2126e
# ╠═099913fe-ac20-4d1c-8755-a7ac4d1da78c
# ╠═e7deeab1-a4c1-46d8-8ac1-2da99e1d18d4
# ╠═dd13511c-54bd-4c66-99b0-2342da4768d1
# ╠═782084b4-6e62-4cad-a81c-cd9d11ab6a1c
# ╟─63442521-9fe1-4319-9f13-82daf14b832f
# ╠═a6502774-b01e-49ae-9e1b-5f58373b8525
# ╠═3b46e566-c7f2-4892-91b1-59bfd91c91b1
# ╠═d6189855-bed3-4758-b14d-b93849cb04b9
# ╟─710d8bdc-32e8-484a-a698-8ff30c1ad27a
# ╠═02529c4d-79d6-4210-8add-ec2375ac2715
# ╠═05fb0941-85ad-4768-9adb-119f15542c1f
# ╠═38573cec-0eb5-47eb-b15f-0dc759f6872f
# ╠═3500d280-3fe2-4054-aed5-dd98491216d5
# ╟─49f63b21-88e5-49b6-b5fc-2708553a9637
# ╟─87ec64a5-ff9c-4219-87e3-dfeddb628684
# ╠═ca2773ba-c781-4786-be3d-d9c0d549a1c8
# ╠═152ac4eb-7eb5-41f8-b9d7-bffa5e61f4c2
# ╠═2016bd95-739a-4227-a4a0-bab46676296c
# ╠═92d7b42d-3ff9-4516-b695-11c873d429ce
# ╠═4650b3c9-c490-4520-a45d-4666b5f31bc8
# ╟─df143390-dcd2-4709-b029-a56fc1e0f73c
# ╠═02aa1aa3-ef47-47be-bbbd-bae4feb73b4b
# ╠═0770f489-ab7b-4fa3-aeb5-295cff0371b4
# ╠═4b1009c0-5f5a-4373-aaed-0b1a1c9a7bc9
# ╠═09f3d176-9529-4fdc-bc74-5578b23dc06f
# ╠═250df9a1-a935-4325-b3f1-448a172fec1a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
