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
	using CSV, DataFrames, CairoMakie, LaTeXStrings
end

# ╔═╡ 2ccdb9d0-9c64-438b-97c2-9684de47efa8
using Statistics, HypothesisTests, MultipleTesting

# ╔═╡ 5aacf9bf-2c83-48b6-94cf-eb6f14b1a163
begin
	using PlutoUI

	_threshold_logFC = @bind threshold_logFC PlutoUI.Slider(0:0.1:3, default = 1.0, show_value = true)
	_threshold_p = @bind threshold_p PlutoUI.Slider(0.0:0.01:0.05, default = 1.0, show_value = true)
	_threshold_FDR  = @bind threshold_FDR PlutoUI.Slider(0.0:0.01:0.5, default = 0.05, show_value = true)

	md"""
	Threshold logFC: $(_threshold_logFC)

	Threshold P $(_threshold_p)
	
	Threshold FDR: $(_threshold_FDR)
	"""
end

# ╔═╡ 7db09fd1-1bd8-4586-a108-e567d6c5ef9b
include("../../../modules/Parsers.jl")


# ╔═╡ f24dd5a5-9adc-412e-9918-3146cd5614ff
md"""
# Space Radiation Transcriptomics
## Importing/Inspecting Data
"""

# ╔═╡ 72c8069f-965f-4517-9d83-bad3b264080c
begin
	availableIDs = readdir("../data/GEO")[2:end]
	_geoID = 	
	md"""
	Select GEO ID from locally available series:
	$(@bind geoID Select(availableIDs))
	"""
end

# ╔═╡ 4c8a58bf-db7a-4259-a442-2755b26d1bf9
md"""
**Working Directory**: $(pwd())
"""

# ╔═╡ d8e61a20-4a9d-4fa0-b8a5-a6f4f9d8a7eb
begin

	geo = "../data/GEO"
	path = joinpath(geo, geoID)

	metadata
	meta = CSV.read(joinpath(path, "metadata_processed.csv"), DataFrame)

	# expression data
	files = filter(name -> endswith(name, "_expression.csv"), readdir(path))
	expressions = [joinpath(path, file) for file in files]

	data = construct_expression_matrix(expressions, meta)

	# container = load_GSE(geoID, basepath = "../data/GEO", parser = parse_GSE64375)
	# data = container.expression
	# meta = container.meta
	# meta.radiation_dose = parse.(Float64, first.(split.(meta.radiation_dose, " ")))

	nothing

end

# ╔═╡ 6de783ff-f4e4-41af-8674-36fb0aacb422
md"""
### Expression Data for $(geoID)
"""

# ╔═╡ 5cc89826-1495-4bb8-9bec-008053238c08
data

# ╔═╡ 9ce994e0-a2d1-4957-86b2-96934b771bc9
md"""
### Metadata for $(geoID)
"""

# ╔═╡ ac524e16-a680-4ee4-b406-84e68838543e
meta

# ╔═╡ c73c494a-a56a-4f52-aab1-3f65f3920856
md"""
## Differential Expression Analysis
"""

# ╔═╡ 1133c36a-6736-40fc-82d6-34c4dd0af994
md"""
### Defining gene list, expression matrix and treatment groups
"""

# ╔═╡ 1a979803-c755-43ab-a692-66cf575ed1f0
begin
	genes = data.ID_REF
	X = Matrix(select(data, Not(:ID_REF)))
	samples = names(select(data, Not(:ID_REF)))
	groups = [meta.group[meta.GSM_ID .== id][1] for id in samples]
	println(genes[1:5])
	println(samples[1:5])
	println(groups[1:5])
end

# ╔═╡ 7c3a5bdb-857c-4ba2-8abd-ac8c19bf562a
begin
	control_idx = findall(==("control"), groups)
	irradiated_idx = findall(!=("control"), groups)
	low_idx = findall(==("low"), groups)
	med_idx = findall(==("med"), groups)
	high_idx = findall(==("high"), groups)

	X_control = X[:, control_idx]
	X_irradiated = X[:, irradiated_idx]
	X_low = X[:, low_idx]
	X_med = X[:, med_idx]
	X_high = X[:, high_idx]
end

# ╔═╡ fba30a5a-8aae-41a9-a376-a79d130b8ec0
begin
	println("Control samples:  ", size(X_control, 2))
	println("Low-dose samples: ", size(X_low, 2))
	println("Med-dose samples: ", size(X_med, 2))
	println("High-dose samples:", size(X_high, 2))
end

# ╔═╡ 83764eb8-cd08-4acb-8658-96e180e36786
md"""
### Quantify Differential Expression


We'll now iterate over each row of the expression matrix $X$ and compute the:
- **$\log_{2}$ fold change** between groups
- **p-value** using an unpaired two-sample **t-test**
"""

# ╔═╡ 9374095c-40c0-4d0e-9fa8-a5b45726d9e6
md"""
#### Checking Sample Variance Uniformity

Before performing differential expression testing, we check whether the variance of gene expression across samples is approximately equal between groups.

This matters because different versions of the t-test make different assumptions:

**Student’s t-test** assumes equal variance and uses the test statistic  
$t = \frac{\bar{x}_1 - \bar{x}_2}{s_p \sqrt{\frac{1}{n_1} + \frac{1}{n_2}}}$  
with pooled variance  
$s_p^2 = \frac{(n_1 - 1)s_1^2 + (n_2 - 1)s_2^2}{n_1 + n_2 - 2}$

**Welch’s t-test** does not assume equal variance:  
$t = \frac{\bar{x}_1 - \bar{x}_2}{\sqrt{ \frac{s_1^2}{n_1} + \frac{s_2^2}{n_2} }}$  
with degrees of freedom:  
$\nu = \frac{(s_1^2/n_1 + s_2^2/n_2)^2}{\frac{(s_1^2/n_1)^2}{n_1 - 1} + \frac{(s_2^2/n_2)^2}{n_2 - 1}}$

By inspecting gene-wise variance between groups, we choose the appropriate test — reducing false positives and improving reliability.
"""


# ╔═╡ 1f3f8f51-f0ed-45a0-ac91-ad34ae6480b4
function check_variance_homoegeneity(X_a, X_b, title::String)
	var_a = mapslices(var, X_a; dims = 2)
	var_b = mapslices(var, X_b; dims = 2)
	ratio = var_a./var_b
	println(title)
	println("Median variance ratio: ", median(ratio))
    println("Max variance ratio: ", maximum(ratio), "\n")
end

# ╔═╡ 08dd95bf-428f-4d97-8042-4ad91dc6671c
begin
	check_variance_homoegeneity(X_control, X_low, "control vs low")
	check_variance_homoegeneity(X_low, X_med, "low vs med")
	check_variance_homoegeneity(X_med, X_high,"med vs high")
	check_variance_homoegeneity(X_high, X_low,"high vs low")
end

# ╔═╡ 9423cb26-99f3-4ea1-b459-964b4db09906
md"""
### Computing DEG statistics

Computes differential expression stats between two groups of samples.

**Inputs**  
- `X_a`, `X_b`: Matrices of gene expression (genes × samples).  
- `method`: Multiple testing correction (default: Benjamini-Hochberg).

**Returns**  
Tuple of vectors `(logFC, pvals, adj_pvals)`:
- `logFC`: log₂(mean(X_a) / mean(X_b)) per gene.  
- `pvals`: Welch’s t-test p-values per gene.  
- `adj_pvals`: Adjusted p-values via `adjust`.

**Logic**  
1. Compute log₂ fold change between group means.  
2. For each gene: run Welch’s t-test unless values are NaN or flat.  
3. Adjust p-values for FDR using `method`.

Use in DEG analysis pipelines for volcano plots, filtering, and enrichment.

"""

# ╔═╡ 65469a91-f2eb-437d-8785-92799bc1ec96
function compute_DEG_stats(X_a::Matrix, X_b::Matrix; method=BenjaminiHochberg())
    logFC = log2.(mean(X_a; dims=2) ./ mean(X_b; dims=2))[:]
    pvals = Float64[]
    for i in 1:size(X_a, 1)
        a, b = X_a[i, :], X_b[i, :]
        if any(isnan, a) || any(isnan, b) || std(a) == 0 || std(b) == 0 || length(a) ≤ 1 || length(b) ≤ 1
            push!(pvals, NaN)
        else
            push!(pvals, pvalue(UnequalVarianceTTest(a, b)))
        end
    end
    adj_pvals = adjust(pvals, method)
    return (logFC, pvals, adj_pvals)
end

# ╔═╡ 6832af54-1711-4be3-a9f0-967d0a63d3e0
begin
	logFC_low, pval_low, adj_low = compute_DEG_stats(X_low, X_control)
	logFC_med, pval_med, adj_med = compute_DEG_stats(X_med, X_control)
	logFC_high, pval_high, adj_high = compute_DEG_stats(X_high, X_control)
end

# ╔═╡ 2c7292da-50c0-4c87-9518-afb13274ff1f
md"""
### Volcano Plot for $(geoID)
"""

# ╔═╡ 420c29de-2c1f-4fbc-9a15-5723e2cb6335
function make_volcano_df(logFC, pval, adj_pval; label)
	DataFrame(
		gene = genes,
		logFC = logFC,
		pval = pval,
		adj_pval = adj_pval,
		neglog10adj = -log10.(adj_pval),
		significant_p = (abs.(logFC) .> threshold_logFC) .&& (pval .< threshold_p),
		significant_adj = (abs.(logFC) .> threshold_logFC) .&& (adj_pval .< threshold_FDR),
		label = label
	)
end

# ╔═╡ 56390a70-f2db-462e-ad6c-e986773f8416
begin
	df_low  = make_volcano_df(logFC_low,  pval_low,  adj_low,  label="low")
	df_med  = make_volcano_df(logFC_med,  pval_med,  adj_med,  label="med")
	df_high = make_volcano_df(logFC_high, pval_high, adj_high, label="high")
	nothing
end

# ╔═╡ bbe8bc0b-af61-4c2a-8911-0699a91c16a5
begin
	df_all = vcat(df_low, df_med, df_high)

	# Pivot to wide format to compare DEGs
	using DataFramesMeta
	
	df_wide = unstack(df_all[df_all.significant_p, [:gene, :label]], :label, :gene)
	
	# Example: genes significant in all doses
	common_genes = intersect(intersect(df_low.gene[df_low.significant_p], df_med.gene[df_med.significant_p]), df_high.gene[df_high.significant_p])

	describe(df_all)

end

# ╔═╡ 06908110-4e5a-4b2e-8670-972e22425182
begin
	dosedict = Dict(0.0 => "control", 0.15 => "low", 0.3 => "med", 1.5 => "high")
	meta.group = get.(Ref(dosedict), meta[!, "radiation_dose"], "unknown")
	DataFrames.combine(groupby(meta, :group), nrow => :count)
end

# ╔═╡ 8c941511-e6ae-4722-9fe7-55eddf709807
md"""
### Visualising DEGs
"""

# ╔═╡ 7c2c78bc-19dc-4630-9b78-f1007f506108
function plot_volcano_grid(df_low, df_med, df_high)
    fig = Figure(size = (1800, 600))  # wider layout

    # Axis 1: Low dose
    ax1 = Axis(fig[1, 1],
        title = "Low Dose",
        xlabel = L"\log_{2}(\mathrm{Fold~Change})",
        ylabel = L"-\log_{10}(\mathrm{UnAdjusted~p\text{-}value})",
        titlesize = 30, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 16, yticklabelsize = 16,
    )

    # Axis 2: Medium dose
    ax2 = Axis(fig[1, 2],
        title = "Medium Dose",
        xlabel = L"\log_{2}(\mathrm{Fold~Change})",
        titlesize = 30, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 16, yticklabelsize = 16,
    )

    # Axis 3: High dose
    ax3 = Axis(fig[1, 3],
        title = "High Dose",
        xlabel = L"\log_{2}(\mathrm{Fold~Change})",
        titlesize = 30, xlabelsize = 22, ylabelsize = 22,
        xticklabelsize = 16, yticklabelsize = 16,
    )

    function safe_scatter!(ax, df)
        x, y, sig = df.logFC, -log10.(df.pval), df.significant_p
        mask = (x, y, s) -> isfinite(x) && isfinite(y) && !ismissing(s)
        all_idx = [i for i in eachindex(x) if mask(x[i], y[i], sig[i])]
        sig_idx = [i for i in all_idx if sig[i]]
        scatter!(ax, x[all_idx], y[all_idx]; color = (:gray, 0.3), markersize = 18)
        scatter!(ax, x[sig_idx], y[sig_idx]; color = :red, markersize = 18)
    end

    safe_scatter!(ax1, df_low)
    safe_scatter!(ax2, df_med)
    safe_scatter!(ax3, df_high)

    fig
end


# ╔═╡ ea3a08c9-21ec-4765-be2b-494803b5b7e2
fig = plot_volcano_grid(df_low, df_med, df_high)

# ╔═╡ 50e0adb4-2312-4da0-b9e3-00d13cc9b43a
save(joinpath("../figures", "$geoID volcano.png"), fig)

# ╔═╡ 21f7fa53-e6c1-4d4a-803a-a3835967a620
df_all[df_all.significant_p, :]


# ╔═╡ Cell order:
# ╠═7b8d09e8-5e1a-489f-bdcc-6900d98f800b
# ╠═7db09fd1-1bd8-4586-a108-e567d6c5ef9b
# ╟─f24dd5a5-9adc-412e-9918-3146cd5614ff
# ╟─72c8069f-965f-4517-9d83-bad3b264080c
# ╟─4c8a58bf-db7a-4259-a442-2755b26d1bf9
# ╟─d8e61a20-4a9d-4fa0-b8a5-a6f4f9d8a7eb
# ╟─6de783ff-f4e4-41af-8674-36fb0aacb422
# ╟─5cc89826-1495-4bb8-9bec-008053238c08
# ╟─9ce994e0-a2d1-4957-86b2-96934b771bc9
# ╟─ac524e16-a680-4ee4-b406-84e68838543e
# ╠═06908110-4e5a-4b2e-8670-972e22425182
# ╟─c73c494a-a56a-4f52-aab1-3f65f3920856
# ╟─1133c36a-6736-40fc-82d6-34c4dd0af994
# ╟─1a979803-c755-43ab-a692-66cf575ed1f0
# ╟─7c3a5bdb-857c-4ba2-8abd-ac8c19bf562a
# ╟─fba30a5a-8aae-41a9-a376-a79d130b8ec0
# ╟─83764eb8-cd08-4acb-8658-96e180e36786
# ╠═2ccdb9d0-9c64-438b-97c2-9684de47efa8
# ╟─9374095c-40c0-4d0e-9fa8-a5b45726d9e6
# ╟─1f3f8f51-f0ed-45a0-ac91-ad34ae6480b4
# ╟─08dd95bf-428f-4d97-8042-4ad91dc6671c
# ╟─9423cb26-99f3-4ea1-b459-964b4db09906
# ╟─65469a91-f2eb-437d-8785-92799bc1ec96
# ╠═6832af54-1711-4be3-a9f0-967d0a63d3e0
# ╟─2c7292da-50c0-4c87-9518-afb13274ff1f
# ╟─5aacf9bf-2c83-48b6-94cf-eb6f14b1a163
# ╟─420c29de-2c1f-4fbc-9a15-5723e2cb6335
# ╠═56390a70-f2db-462e-ad6c-e986773f8416
# ╠═bbe8bc0b-af61-4c2a-8911-0699a91c16a5
# ╟─8c941511-e6ae-4722-9fe7-55eddf709807
# ╟─7c2c78bc-19dc-4630-9b78-f1007f506108
# ╠═ea3a08c9-21ec-4765-be2b-494803b5b7e2
# ╠═50e0adb4-2312-4da0-b9e3-00d13cc9b43a
# ╠═21f7fa53-e6c1-4d4a-803a-a3835967a620
