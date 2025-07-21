### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ 7b8d09e8-5e1a-489f-bdcc-6900d98f800b
begin
	using Pkg
	Pkg.activate("/Users/fomo/Documents/Kaizen/code/bioinformatics")
	Pkg.status()
	using CSV, DataFrames
end

# ╔═╡ 2ccdb9d0-9c64-438b-97c2-9684de47efa8
using Statistics, HypothesisTests, MultipleTesting

# ╔═╡ 5aacf9bf-2c83-48b6-94cf-eb6f14b1a163
begin
	using PlutoUI

	_threshold_logFC = @bind threshold_logFC Slider(0:0.1:3, default = 1.0, show_value = true)
	_threshold_FDR  = @bind threshold_FDR Slider(0.0:0.01:0.5, default = 0.05, show_value = true)

	md"""
	Threshold logFC: $(_threshold_logFC)
	
	Threshold FDR: $(_threshold_FDR)
	"""
end

# ╔═╡ 7db09fd1-1bd8-4586-a108-e567d6c5ef9b
begin
	include("../../../modules/GeoPreprocess.jl")
	nothing
end

# ╔═╡ f24dd5a5-9adc-412e-9918-3146cd5614ff
md"""
# Space Radiation Transcriptomics
## Importing/Inspecting Data
"""

# ╔═╡ 39ed2cf9-95e7-4803-89c8-18d9565a4c87


# ╔═╡ d8e61a20-4a9d-4fa0-b8a5-a6f4f9d8a7eb
begin
	
	geo = "projects/space/data/GEO/GSE63952"

	# metadata
	meta = CSV.read(joinpath(geo, "metadata_processed.csv"), DataFrame)

	# expression data
	files = filter(name -> endswith(name, "_expression.csv"), readdir(geo))
	expressions = [joinpath(geo, file) for file in files]

	data = construct_expression_matrix(expressions, meta)

	first(data, 5), first(meta, 5)

end

# ╔═╡ 1133c36a-6736-40fc-82d6-34c4dd0af994
md"""
#### Defining gene list, expression matrix and treatment groups
"""

# ╔═╡ 1a979803-c755-43ab-a692-66cf575ed1f0
begin
	genes = data.ID_REF
	X = Matrix(select(data, Not(:ID_REF)))
	samples = names(select(data, Not(:ID_REF)))
	groups = [meta.group[meta.GSM_ID .== id][1] for id in samples] 
end

# ╔═╡ c73c494a-a56a-4f52-aab1-3f65f3920856
md"""
## Differential Expression Analysis
"""

# ╔═╡ 7c3a5bdb-857c-4ba2-8abd-ac8c19bf562a
begin
	control_idx = findall(==("control"), groups)
	irradiated_idx = findall(==("irradiated"), groups)

	X_control = X[:, control_idx]
	X_irradiated = X[:, irradiated_idx]
end

# ╔═╡ 83764eb8-cd08-4acb-8658-96e180e36786
md"""
### Quantify Differential Expression


We'll now iterate over each row of the expression matrix $X$ and compute the:
- **$\log_{2}$ fold change** between groups
- **p-value** using an unpaired two-sample **t-test**
"""

# ╔═╡ 65469a91-f2eb-437d-8785-92799bc1ec96
logFC = log2.(mean(X_irradiated; dims = 2) ./ mean(X_control; dims = 2))[:]

# ╔═╡ 7066f880-8ca5-4a84-9b15-7d951b01cfec
pval = [pvalue(UnequalVarianceTTest(X_irradiated[i, :], X_control[i, :])) for i in 1:size(X, 1)]

# ╔═╡ e50bfc67-606f-4ae3-91ac-1475cd30ef9a
adj_pval_BH = adjust(pval, BenjaminiHochberg()) 

# ╔═╡ 7e878459-3aea-4474-aed0-236d04208c30
adj_pval_BON = adjust(pval, Bonferroni()) 

# ╔═╡ 657f8f19-4c24-418c-b99a-e47c4022af7d


# ╔═╡ c9de8cbf-9160-4568-a525-3d3b3c89cabf
volcano_df = DataFrame(
    gene = genes,
    logFC = logFC,
    adj_pval_BH = adj_pval_BH,
    adj_pval_BON = adj_pval_BON,
    neglog10FDR = -log10.(adj_pval_BH),
    neglog10BON = -log10.(adj_pval_BON),
    significant_BH = (abs.(logFC) .> threshold_logFC) .&& (adj_pval_BH .< threshold_FDR),
    significant_BON = (abs.(logFC) .> threshold_logFC) .&& (adj_pval_BON .< threshold_FDR)
)

# ╔═╡ b18a68c7-a8a0-4171-8017-db22a4b1cf0b
begin

	using CairoMakie
	
	x = volcano_df.logFC
	y = volcano_df.neglog10FDR
	is_BH = volcano_df.significant_BH
	is_BON = volcano_df.significant_BON

	fig = Figure(size = (1500, 1000))
	ax = Axis(
		fig[1,1],

		# TITLE
		title = "Volcano Plot GSE63952",
		titlesize = 49,
        titlegap = 21,

		# X-LABEL
		xlabel = L"\log_{2}(\mathrm{Fold~Change})",
		xlabelsize = 30,
		xticklabelsize = 21,
		xlabelpadding = 21,

		# Y-LABEL
		ylabel = L"-\log_{10}(\mathrm{Adjusted~p\text{-}value})",	
        ylabelsize = 30,
        yticklabelsize = 21,
		ylabelpadding = 21,

		# SPINES
        spinewidth = 1.5,
        leftspinevisible = true,
        bottomspinevisible = true,
        rightspinevisible = false,
        topspinevisible = false,

		# GRIDS
        xgridvisible = true,
		xminorgridvisible = true,
        ygridvisible = true,
		yminorgridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xgridwidth = 2.1,
		ygridwidth = 2.1,
        xgridcolor = RGBAf(0.7, 0.7, 0.7, 0.7),
		ygridcolor = RGBAf(0.7, 0.7, 0.7, 0.7)
		
	)
	scatter!(ax, x, y, color = (:gray, 0.3), markersize = 21, )
	scatter!(ax, x[is_BH], y[is_BH], color = :red, markersize = 21)
	scatter!(ax, x[is_BON], y[is_BON], color = :purple, markersize = 21)

end

# ╔═╡ 5019214b-44de-4880-8977-7d6b868060c6
volcano_df[(volcano_df.significant_BH .| volcano_df.significant_BON), :]

# ╔═╡ 8c941511-e6ae-4722-9fe7-55eddf709807
md"""
### Visualising DEGs
"""

# ╔═╡ 9d33a328-3cef-4026-8302-dfb1be2e4289
fig

# ╔═╡ Cell order:
# ╠═7b8d09e8-5e1a-489f-bdcc-6900d98f800b
# ╠═7db09fd1-1bd8-4586-a108-e567d6c5ef9b
# ╟─f24dd5a5-9adc-412e-9918-3146cd5614ff
# ╠═39ed2cf9-95e7-4803-89c8-18d9565a4c87
# ╠═d8e61a20-4a9d-4fa0-b8a5-a6f4f9d8a7eb
# ╠═1133c36a-6736-40fc-82d6-34c4dd0af994
# ╠═1a979803-c755-43ab-a692-66cf575ed1f0
# ╟─c73c494a-a56a-4f52-aab1-3f65f3920856
# ╠═7c3a5bdb-857c-4ba2-8abd-ac8c19bf562a
# ╠═83764eb8-cd08-4acb-8658-96e180e36786
# ╠═2ccdb9d0-9c64-438b-97c2-9684de47efa8
# ╠═65469a91-f2eb-437d-8785-92799bc1ec96
# ╠═7066f880-8ca5-4a84-9b15-7d951b01cfec
# ╠═e50bfc67-606f-4ae3-91ac-1475cd30ef9a
# ╠═7e878459-3aea-4474-aed0-236d04208c30
# ╠═5aacf9bf-2c83-48b6-94cf-eb6f14b1a163
# ╠═657f8f19-4c24-418c-b99a-e47c4022af7d
# ╠═c9de8cbf-9160-4568-a525-3d3b3c89cabf
# ╠═5019214b-44de-4880-8977-7d6b868060c6
# ╟─8c941511-e6ae-4722-9fe7-55eddf709807
# ╠═b18a68c7-a8a0-4171-8017-db22a4b1cf0b
# ╠═9d33a328-3cef-4026-8302-dfb1be2e4289
