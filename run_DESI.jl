using FITSIO, Statistics, DataFrames, JSON, DataStructures, CSV, Revise
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer
using LinesInTheSky


function analyze_single_spec(input_path, output_path, row)
    mkpath("$(output_path)/JSON")
    mkpath("$(output_path)/HTML")

    input_filename = "$(input_path)/TXT/$(row[:id_DESI_DR1])_$(row[:Z]).txt"
    output_filename = "$(output_path)/JSON/$(row[:id_DESI_DR1])_$(row[:Z]).json"

    if !isfile(input_filename)
        error("Input file $input_filename do not exists.")
    end
    if isfile(output_filename)
        println("Output file $output_filename already exists. Skipping.")
        return
    end

    println(); println()
    @info "Analyzing file $(input_filename)"
    spec = Spectrum(Val(:ASCII), input_filename, columns=[1,2,3], resolution= 2857, label=string(row[:id_DESI_DR1]))
    recipe = CRecipe{WP9Type1IR}(redshift=row[:Z], use_host_template=false, Av=0.0)
    res = analyze(recipe, spec)
    GModelFit.serialize(output_filename, res.bestfit, res.fsumm)
    GModelFitViewer.serialize_html(filename="$(output_path)/HTML/$(row[:id_DESI_DR1])_$(row[:Z]).html", res)

    # Write additional info in a JSON file
    aux = Dict{Symbol, Any}()
    aux[:SNR] = median(abs.(values(res.data) ./ uncerts(res.data)))
    f = open("$(output_path)/JSON/$(row[:id_DESI_DR1])_$(row[:Z])_aux.json", "w")
    write(f, JSON.json(aux))
    close(f)

    return res
end


function run_analysis(input_path, output_path, catalog)
    Threads.@threads for i in 1:nrow(catalog)
		try
		    analyze_single_spec(input_path, output_path, catalog[i, :])
		catch err
		    display(err)
		    println("Failed to fit spectrum for $(catalog[i, :id_DESI_DR1]).")
		end
    end
end


function read_results(output_path, catalog)
    out = DataFrame(ID_DESI=Int[], ID_EUCLID=Int[], Redshift=Float64[], redchisq=Float64[], SNR=Float64[], NPOINTS=Float64[], Html_serial=String[])
    for i in 1:nrow(catalog)
        filename = "$(output_path)/JSON/$(catalog[i, :id_DESI_DR1])_$(catalog[i,:Z]).json"
        if isfile(filename)
            @info "Reading $filename ..."
            bestfit, fsumm = GModelFit.deserialize(filename)
            aux = JSON.Parser.parsefile("$(output_path)/JSON/$(catalog[i, :id_DESI_DR1])_$(catalog[i, :Z])_aux.json")
            push!(out, [catalog[i, :id_DESI_DR1], catalog[i, :object_id], catalog[i, :Z], fsumm.fitstat, aux["SNR"],fsumm.ndata, "", fill(missing, ncol(out)-7)...])
            out[end, :Html_serial] = "$(output_path)/HTML/$(catalog[i, :id_DESI_DR1])_$(catalog[i, :Z]).html"
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


input_path  =   "input_DESI"
output_path = "results_DESI"

# Read input catalog
f = FITS("$(input_path)/Q1_DESI_DR1_QSOAGNVAC_QSO.fits")
catalog = DataFrame(f[2])
close(f)

# Run analysis
run_analysis(input_path, output_path, catalog)

# Read results from JSON files
results = read_results(output_path, catalog)

# Write results in a FITS file
data = OrderedDict{String, Vector}()
for cname in names(results)
    col = results[:, cname]
    typ = nonmissingtype(eltype(col))
    data[cname] = replace(col, missing => (typ <: Number  ?  -1  :  ""))
end
f = FITS("$(output_path)/QSFIT_RESULTS.fits", "w")
write(f, data)
close(f)
