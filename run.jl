using FITSIO, Statistics, DataFrames, JSON, DataStructures, CSV, Revise
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer
using LinesInTheSky

function analyze_single_spec(row)
    mkpath("results/JSON")
    mkpath("results/HTML")

    input_filename = "input/spectra/$(row[:object_id]).txt"
    output_filename = "results/JSON/$(row[:object_id]).json"

    if !isfile(input_filename)
        error("Input file $input_filename do not exists.")
    end
    if isfile(output_filename)
        println("Output file $output_filename already exists. Skipping.")
        return
    end

    println(); println()
    @info "Analyzing file $(input_filename)"
    spec = Spectrum(Val(:ASCII), input_filename, columns=[1,2,5], label=string(row[:object_id]), resolution = 450)
    recipe = CRecipe{WP9Type1IR}(redshift=row[:Z], use_host_template=true, Av=0.0, n_nuisance=2)
    res = analyze(recipe, spec)
    GModelFit.serialize(output_filename, res.bestfit, res.fsumm)
    GModelFitViewer.serialize_html(filename="results/HTML/$(row[:object_id]).html", res)

    # Write additional info in a JSON file
    aux = Dict{Symbol, Any}()
    aux[:SNR] = median(abs.(values(res.data) ./ uncerts(res.data)))
    f = open("results/JSON/$(row[:object_id])_aux.json", "w")
    write(f, JSON.json(aux))
    close(f)
end


function run_analysis(catalog)
    Threads.@threads for i in 1:nrow(catalog)
        try
            analyze_single_spec(catalog[i, :])
        catch err
            display(err)
            println("Failed to fit spectrum for $(catalog[i, :object_id]).")
        end
    end
end



function read_results(catalog)
    out = DataFrame(ID=Int[], Redshift=Float64[],Source=String[], redchisq=Float64[], SNR=Float64[], NPOINTS=Float64[], Html_serial=String[])
    for i in 1:nrow(catalog)
        filename = "results/JSON/$(catalog[i, :object_id]).json"
        if isfile(filename)
            @info "Reading $filename ..."
            bestfit, fsumm = GModelFit.deserialize(filename)
            aux = JSON.Parser.parsefile("results/JSON/$(catalog[i, :object_id])_aux.json")
            push!(out, [catalog[i, :object_id], catalog[i, :Z], catalog[i, :CAT], fsumm.fitstat, aux["SNR"],fsumm.ndata, "", fill(missing, ncol(out)-7)...])
            out[end, :Html_serial] = "results/HTML/$(catalog[i, :object_id]).html"
            for (cname, comp) in bestfit
                 for (pname, par) in comp
                     colname = Symbol(cname, :_, pname)
                     if !(string(colname) in names(out))
                         out[!,        colname        ] = missings(Float64, nrow(out))
                         out[!, Symbol(colname, :_unc)] = missings(Float64, nrow(out))
                     end

                     if isnothing(par.patch)
                         out[end,        colname        ] = par.val
                         out[end, Symbol(colname, :_unc)] = par.unc
                     else
                         out[end,        colname        ] = par.actual
                         out[end, Symbol(colname, :_unc)] = NaN
                     end
                 end
            end
        end
    end

    # Set Html_serial as last column
    select!(out, [filter(x -> x != "Html_serial", names(out)); "Html_serial"])
    return out
end

# Read input catalog
f = FITS("input/catalog.fits")
catalog = DataFrame(f[2])
close(f)

# Run analysis
run_analysis(catalog)


# Read results from JSON files
results = read_results(catalog)

# Write results in a FITS file
data = OrderedDict{String, Vector}()
for cname in names(results)
    col = results[:, cname]
    typ = nonmissingtype(eltype(col))
    data[cname] = replace(col, missing => (typ <: Number  ?  -1  :  ""))
end
f = FITS("results/QSFIT_RESULTS.fits", "w")
write(f, data)
close(f)
