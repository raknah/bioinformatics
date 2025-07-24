function parse_GSE64375(df::DataFrame)
    # Prepare a Vector of Dicts, one per row
    rows = Vector{Dict{Symbol,Any}}(undef, nrow(df))
    
    for (i, row) in enumerate(eachrow(df))
        # Start each record with the GSM_ID
        d = Dict{Symbol,Any}(:GSM_ID => row.GSM_ID)
        
        # Split the “Characteristics” string on ‘;’
        for part in split(row.Characteristics, ";")
            # Split each part into key and value (at most 2 pieces)
            kv = split(part, ":", limit = 2)
            if length(kv) == 2
                # Normalize the key to a Symbol (spaces → underscores)
                key = Symbol(replace(strip(kv[1]), r"\s+" => "_"))
                # Strip whitespace from the value
                d[key] = strip(kv[2])
            end
        end
        
        rows[i] = d
    end
    
    # Convert the array of Dicts into a DataFrame
    return DataFrame(rows)
end