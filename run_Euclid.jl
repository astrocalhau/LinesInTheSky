using FITSIO, Statistics, DataFrames, JSON, DataStructures, CSV, Revise, Dierckx
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer
using LinesInTheSky

include("common_functs.jl")

function analyze_single_spec(input_path, output_path, row)
    mkpath("$(output_path)/JSON")
    mkpath("$(output_path)/HTML")

    input_filename = "$(input_path)/spectra/$(row[:object_id]).txt"
    output_filename = "$(output_path)/JSON/$(row[:object_id]).json"

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
    recipe = CRecipe{WP9Type1IR}(redshift=row[:Z], use_host_template=false, Av=0.0, n_nuisance=2)
    resNoHost = analyze(recipe, spec)

    recipe.use_host_template = true
    resWithHost = analyze(recipe, spec)

    if resNoHost.fsumm.fitstat < resWithHost.fsumm.fitstat
        res = resNoHost
    else
        res = resWithHost
    end

    GModelFit.serialize(output_filename, res.bestfit, res.fsumm, res.data)
    GModelFitViewer.serialize_html(filename="$(output_path)/HTML/$(row[:object_id]).html", res)

    # Write additional info in a JSON file
    aux = Dict{Symbol, Any}()
    aux[:SNR] = median(abs.(values(res.data) ./ uncerts(res.data)))
    aux[:L3000] = ((res.post[:Quality_flags][:QSOcont] == 0)  ?  cont_lambdaLlambda(res.bestfit, 3000.)  :  NaN)
    aux[:L5100] = ((res.post[:Quality_flags][:QSOcont] == 0)  ?  cont_lambdaLlambda(res.bestfit, 5100.)  :  NaN)
    aux[:reliable] = Symbol[]
    for cname in keys(res.bestfit)
        if res.post[:Quality_flags][cname] == 0
            push!(aux[:reliable], cname)
        end
    end
    f = open("$(output_path)/JSON/$(row[:object_id])_aux.json", "w")
    write(f, JSON.json(aux, allownan=true))
    close(f)

    return res
end


function run_analysis(input_path, output_path, catalog)
    Threads.@threads for i in 1:nrow(catalog)
        try
            analyze_single_spec(input_path, output_path, catalog[i, :])
        catch err
            display(err)
            println("Failed to fit spectrum for $(catalog[i, :object_id]).")
        end
    end
end


function read_results(output_path, catalog)
    out = DataFrame(ID=Int[], Redshift=Float64[],Source=String[], redchisq=Float64[], NPOINTS=Float64[],
                    SNR=Float64[], L3000=Float64[], L5100=Float64[], Html_serial=String[])
    allowmissing!(out, [:L3000, :L5100])
    ncol_initial = ncol(out)
    for i in 1:nrow(catalog)
        filename = "$(output_path)/JSON/$(catalog[i, :object_id]).json"
        if isfile(filename)
            @info "Reading $filename ..."
            bestfit, fsumm = GModelFit.deserialize(filename)
            aux = JSON.parsefile("$(output_path)/JSON/$(catalog[i, :object_id])_aux.json", allownan=true)
            push!(out, [catalog[i, :object_id], catalog[i, :Z], catalog[i, :CAT], fsumm.fitstat, fsumm.ndata,
                        aux["SNR"], aux["L3000"], aux["L5100"], "", fill(missing, ncol(out)-ncol_initial)...])
            out[end, :Html_serial] = "$(output_path)/HTML/$(catalog[i, :object_id]).html"
            for (cname, comp) in bestfit
                (cname in [:QSOcont, :Galaxy, :Ironuv, :Ironoptbr, :Ironoptna,
                           :Ha_br, :Ha_na, :Hb_br, :Hb_na, :Pab_br, :HeI_10832_br,
                           :MgII_2798_br, :OIII_4959, :OIII_5007, :OIII_5007_bw])  ||  continue

                if !(string(cname) * "_reliable" in names(out))
                    out[!, Symbol(cname, :_reliable)] = missings(Int64, nrow(out))
                end
                out[end, Symbol(cname, :_reliable)] = ((string(cname) in aux["reliable"])  ?  1  :  0)
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


input_path  =   "input/input_Euclid"
output_path = "results_Euclid"

# Read input catalog
f = FITS("$(input_path)/catalog.fits")
catalog = DataFrame(f[2])
close(f)

# Run analysis
run_analysis(input_path, output_path, catalog)

# Read results from JSON files
results = read_results(output_path, catalog)


# Calculates Mbh
add_MBH_Hb_WuShen2022!(results)
add_MBH_MgII_WuShen2022!(results)
add_MBH_Ha_ShenLiu2012!(results)
add_MBH_Ha_Ricci!(results)
add_MBH_Hb_Ricci!(results)
add_MBH_Ha_L5100_Ricci!(results)
add_MBH_MgII_Ricci!(results)
add_MBH_Pab_Ricci!(results)
add_MBH_HeI_Ricci!(results)

results.MBH_mean .= NaN
for i in 1:nrow(results)
    results[i, :MBH_mean] = mean_handle_NaN([results[i, :MBH_Hb_WuShen2022], results[i, :MBH_MgII_WuShen2022], results[i, :MBH_Ha_ShenLiu2012]])
end

# Calculates Lbol and Eddington ratios
add_Lbol_eddratio!(results)


# Creates Quality cut columns
function snr_min_at_redchisq(x; redchisq=3.5, SNR=3)
    a = log10.([redchisq, SNR])
    lx = log10(x)
    (lx < a[1])  &&  (return 10^a[2])
    lx -= a[1]
    return 10^(0.6 * lx + a[2])
end

results[!, :qcut] = ((results.NPOINTS .> 450)                                        .&
                     (results.SNR .> snr_min_at_redchisq.(results.redchisq))         .&
                     (results.QSOcont_alpha .> -5) .& (results.QSOcont_alpha .<  5)  .&
                     (results.QSOcont_norm .> results.QSOcont_norm_unc))

# Write results in a FITS file
write_fits("$(output_path)/QSFIT_RESULTS.fits", results)
