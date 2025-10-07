using FITSIO, Statistics, DataFrames, JSON, DataStructures, CSV, Revise, Dierckx
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer
using LinesInTheSky


function cont_lambdaLlambda(model, wavelength)
    try
        return Dierckx.Spline1D(coords(domain(model)), model(:QSOcont), k=1, bc="error")(wavelength) * wavelength * 1e-2
    catch
        return NaN
    end
end


function analyze_single_spec(input_path, output_path, row)
    mkpath("$(output_path)/JSON")
    mkpath("$(output_path)/HTML")

    input_filename = "$(input_path)/TXT/$(row[:id_DESI_DR1]).txt"
    output_filename = "$(output_path)/JSON/$(row[:id_DESI_DR1]).json"

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
    recipe = CRecipe{Type1}(redshift=row[:Z], use_host_template=true, Av=0.0)
    res = analyze(recipe, spec)
    GModelFit.serialize(output_filename, res.bestfit, res.fsumm)
    GModelFitViewer.serialize_html(filename="$(output_path)/HTML/$(row[:id_DESI_DR1]).html", res)

    # Write additional info in a JSON file
    aux = Dict{Symbol, Any}()
    aux[:SNR] = median(abs.(values(res.data) ./ uncerts(res.data)))
    aux[:L3000] = cont_lambdaLlambda(res.bestfit, 3000.)
    aux[:L5100] = cont_lambdaLlambda(res.bestfit, 5100.)
    f = open("$(output_path)/JSON/$(row[:id_DESI_DR1])_aux.json", "w")
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
    out = DataFrame(ID_DESI=Int[], ID_EUCLID=Int[], Redshift=Float64[], redchisq=Float64[], NPOINTS=Float64[],
                    SNR=Float64[], L3000=Float64[], L5100=Float64[], Html_serial=String[])
    allowmissing!(out, [:L3000, :L5100])
    ncol_initial = ncol(out)
    for i in 1:nrow(catalog)
        filename = "$(output_path)/JSON/$(catalog[i, :id_DESI_DR1]).json"
        if isfile(filename)
            @info "Reading $filename ..."
            bestfit, fsumm = GModelFit.deserialize(filename)
            aux = JSON.Parser.parsefile("$(output_path)/JSON/$(catalog[i, :id_DESI_DR1])_aux.json")
            for k in ["L3000", "L5100"]
                isnothing(aux[k])  &&  (aux[k] = missing)
            end
            push!(out, [catalog[i, :id_DESI_DR1], catalog[i, :object_id], catalog[i, :Z], fsumm.fitstat, fsumm.ndata,
                        aux["SNR"], aux["L3000"], aux["L5100"], "", fill(missing, ncol(out)-ncol_initial)...])
            out[end, :Html_serial] = "$(output_path)/HTML/$(catalog[i, :id_DESI_DR1]).html"
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


input_path  =   "input/input_DESI"
output_path = "results_DESI"

# Read input catalog
f = FITS("$(input_path)/Q1_DESI_DR1_QSOAGNVAC_QSO.fits")
catalog = DataFrame(f[2])
close(f)

# Run analysis
run_analysis(input_path, output_path, catalog)

# Read results from JSON files
results = read_results(output_path, catalog)


# Calculates Mbh
results.MBH_Hb_WuShen2022   = 0.91 .+ 0.5  .* log10.(results.L5100) .+ 2 .* log10.(results.Hb_br_fwhm)
results.MBH_MgII_WuShen2022 = 0.74 .+ 0.62 .* log10.(results.L3000) .+ 2 .* log10.(results.MgII_2798_br_fwhm)

# Ha_norm = results.Ha_br_norm
# i = findall(.!ismissing.(results.Ha_na_norm))
# Ha_norm[i] .+= results.Ha_na_norm
# results.MBH_Ha_ShenLiu2012 = 2.216 .+ 0.564 .* log10.((Ha_norm) .* 1e-2) .+ 1.821 .* log10.(results.Ha_br_fwhm)
results.MBH_Ha_ShenLiu2012 .= missing

results.MBH_mean .= 0.
allowmissing!(results, :MBH_mean)
for i in 1:nrow(results)
    try
	    results[i, :MBH_mean] = mean(skipmissing([results[i, :MBH_Hb_WuShen2022], results[i, :MBH_MgII_WuShen2022], results[i, :MBH_Ha_ShenLiu2012]]))
    catch
        results[i, :MBH_mean] = missing
    end
end


# Calculates Lbol
results.Lbol_3000 = 5.15e44 .* results.L3000
results.Lbol_5100 = 9.26e44 .* results.L5100
results.Lbol_mean .= 0.
allowmissing!(results, :Lbol_mean)
for i in 1:nrow(results)
    try
	    results[i, :Lbol_mean] = mean(skipmissing([results[i, :Lbol_3000], results[i, :Lbol_5100]]))
    catch
        results[i, :Lbol_mean] = missing
    end        
end


# Calculates Eddington ratios
results.Ledd_mean = 1.26e38 * 10 .^results.MBH_mean
results.Edd_ratio = results.Lbol_mean ./ results.Ledd_mean


# Write results in a FITS file
data = OrderedDict{String, Vector}()
for cname in names(results)
    col = results[:, cname]
    typ = nonmissingtype(eltype(col))
    if typ <: AbstractFloat
        data[cname] = replace(col, missing => NaN)
    elseif typ <: Integer
        data[cname] = replace(col, missing => -1)
    elseif typ <: String
        data[cname] = replace(col, missing => "")
    else
        error("Unsupported data type: $(typ)")
    end
end
f = FITS("$(output_path)/QSFIT_RESULTS.fits", "w")
write(f, data)
close(f)
