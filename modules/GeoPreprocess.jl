# module GeoPreprocess

using CSV, DataFrames

# ————————————————————————————————————————————————————————————
# Data container
# ————————————————————————————————————————————————————————————

Base.@kwdef struct GSEData
    name::String
    meta::DataFrame
    expression::DataFrame
end

# ————————————————————————————————————————————————————————————
# 1. Expression-matrix builder 
# ————————————————————————————————————————————————————————————

function construct_expression_matrix(paths::Vector{String}, metadata::DataFrame)
    path_map = Dict(split(basename(p), "_")[1] => p for p in paths)
    sample_ids = metadata.GSM_ID
    df0 = CSV.read(path_map[sample_ids[1]], DataFrame)
    gene_ids = df0.ID_REF
    n_genes, n_samples = size(df0, 1), length(sample_ids)
    values = zeros(Float64, n_genes, n_samples)
    for (i, gsm) in enumerate(sample_ids)
        df = CSV.read(path_map[gsm], DataFrame)
        values[:, i] = df.VALUE
    end
    X = DataFrame(values, Symbol.(sample_ids))
    X.ID_REF = gene_ids
    return X
end

# ————————————————————————————————————————————————————————————
# 2. Generic metadata preprocessor
# ————————————————————————————————————————————————————————————

"""
    preprocess_metadata(name::String;
                        basepath::String="…/GEO",
                        processed_name::String="metadata_processed.csv",
                        parser::Function)

Reads raw metadata.csv, calls `parser(raw_df)` → a cleaned DataFrame,
caches to `processed_name`, and returns it.
"""
function preprocess_metadata(name::String;
                             basepath::String="projects/space/data/GEO",
                             processed_name::String="metadata_processed.csv",
                             parser::Function)
    path        = joinpath(basepath, name)
    raw_path    = joinpath(path, "metadata.csv")
    out_path    = joinpath(path, processed_name)

    if isfile(out_path)
        return CSV.read(out_path, DataFrame)
    end

    raw = CSV.read(raw_path, DataFrame)
    meta = parser(raw)                              # ← your custom parser
    CSV.write(out_path, meta)
    return meta
end

# ————————————————————————————————————————————————————————————
# 3. Generic loader that ties it all together
# ————————————————————————————————————————————————————————————

"""
    load_GSE(name::String;
             basepath::String="…/GEO",
             parser::Function)

Loads metadata via `preprocess_metadata(..., parser=parser)`,
then builds the expression matrix, and returns `GSEData`.
"""
function load_GSE(name::String;
                  basepath::String="projects/space/data/GEO",
                  parser::Function)
    meta       = preprocess_metadata(name;
                                     basepath=basepath,
                                     parser=parser)
    expr_files = readdir(joinpath(basepath, name)) .|> f -> joinpath(basepath,name,f)
    expr_paths = filter(f->endswith(f, "_expression.csv"), expr_files)
    expr       = construct_expression_matrix(expr_paths, meta)

    return GSEData(name=name, meta=meta, expression=expr)
end

# end # module
