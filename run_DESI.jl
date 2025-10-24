using Revise
using Base.Threads, FITSIO, DataFrames, DataStructures, Printf, Statistics, StatsBase
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer
using LinesInTheSky

include("common_functs.jl")

function analyze_single_spec(input_path, output_path, row; clob=false)
    mkpath("$(output_path)/JSON")
    mkpath("$(output_path)/HTML")

    input_filename = "$(input_path)/TXT/$(row[:id_DESI_DR1]).txt"
    output_filename = "$(output_path)/JSON/$(row[:id_DESI_DR1]).json"

    if !isfile(input_filename)
        error("Input file $input_filename do not exists.")
    end
    if isfile(output_filename)  &&  !clob
        println("Output file $output_filename already exists. Skipping.")
        return
    end

    println(); println()
    @info "Analyzing file $(input_filename)"
    spec = Spectrum(Val(:ASCII), input_filename, columns=[1,2,3], resolution= 2857, label=string(row[:id_DESI_DR1]))
    recipe = CRecipe{Type1}(redshift=row[:Z], use_host_template=true, Av=0.0)
    res = analyze(recipe, spec)
    QSFit.serialize(output_filename, res)
    output_html = replace(output_filename, "JSON" => "HTML", "json" => "html")
    @info output_html
    GModelFitViewer.serialize_html(filename=output_html, res)
    return res
end


function run_analysis(input_path, output_path, catalog)
    # bash -c 'echo "`ls results/HTML | wc -w` / `ls input/ | wc -w`" | bc -l'
    SENTINEL = -1
    channel = Channel{Int64}(30)

    function consumer(channel::Channel)
        while true
            i = take!(channel)
            if i == SENTINEL
                put!(channel, SENTINEL) # tell other threads to quit
                break                   # quit this thread
            end
            try
                analyze_single_spec(input_path, output_path, catalog[i, :])
            catch err
                display(err)
                println("Failed to fit spectrum for $(catalog[i, :id_DESI_DR1]).")
            end
        end
    end

    tasks = [@spawn consumer(channel) for i in 1:nthreads()]
    put!.(Ref(channel), 1:nrow(catalog))  # tell threads which row to analyze
    put!(     channel , SENTINEL)         # tell threads to quit
    wait.(tasks)                          # wait for threads to terminate
end


function read_results(output_path, catalog)
    out = [DataFrame(ID_DESI=Int[], ID_EUCLID=Int[], Redshift=Float64[], redchisq=Float64[], NPOINTS=Float64[],
                    SNR=Float64[], L3000=Float64[], L5100=Float64[], Html_serial=String[]) for i in 1:(Threads.nthreads(:interactive) .+ Threads.nthreads(:default))]
    ncol_initial = ncol(out[1])
    Threads.@threads for i in 1:nrow(catalog)
        df = out[Threads.threadid()]
        filename = "$(output_path)/JSON/$(catalog[i, :id_DESI_DR1]).json"
        if isfile(filename)
            @info "Reading $filename ..."
            res = QSFit.deserialize(filename)
            push!(df, [catalog[i, :id_DESI_DR1], catalog[i, :object_id], catalog[i, :Z], res.fsumm.fitstat, res.fsumm.ndata,
                        res.post[:Data_stats][:SNR],
                        res.post[:Continuum_luminosity][:l3000],
                        res.post[:Continuum_luminosity][:l5100],
                        "", fill(missing, ncol(df)-ncol_initial)...])
            df[end, :Html_serial] = "$(output_path)/HTML/$(catalog[i, :id_DESI_DR1]).html"
            for (cname, comp) in res.bestfit
                (cname in [:QSOcont, :Galaxy, :Ironuv, :Ironoptbr, :Ironoptna,
                           :Ha_br, :Ha_na, :Hb_br, :Hb_na, :Pab_br, :HeI_10832_br,
                           :MgII_2798_br, :OIII_4959, :OIII_5007, :OIII_5007_bw])  ||  continue

                if !(string(cname) * "_reliable" in names(df))
                    df[!, Symbol(cname, :_reliable)] = missings(Int64, nrow(df))
                end
                df[end, Symbol(cname, :_reliable)] = ((string(cname) in keys(res.post[:Issues]))  ?  0  :  1)
                for (pname, par) in comp
                    colname = Symbol(cname, :_, pname)
                    if !(string(colname) in names(df))
                        df[!,        colname        ] = missings(Float64, nrow(df))
                        df[!, Symbol(colname, :_unc)] = missings(Float64, nrow(df))
                    end

                    if isnothing(par.patch)
                        df[end,        colname        ] = par.val
                        df[end, Symbol(colname, :_unc)] = par.unc
                    else
                        df[end,        colname        ] = par.actual
                        df[end, Symbol(colname, :_unc)] = NaN
                    end
                end
            end

            for assoc in [:Ha_br_assoc, :Hb_br_assoc]
                if assoc in keys(res.post)
                    if !("$(assoc)_norm" in names(df))
                        df[!, Symbol(assoc, :_norm)] = missings(Float64, nrow(df))
                        df[!, Symbol(assoc, :_fwhm)] = missings(Float64, nrow(df))
                        df[!, Symbol(assoc, :_voff)] = missings(Float64, nrow(df))
                    end
                    df[end, Symbol(assoc, :_norm)] = res.post[assoc][:norm]
                    df[end, Symbol(assoc, :_fwhm)] = res.post[assoc][:fwhm]
                    df[end, Symbol(assoc, :_voff)] = res.post[assoc][:voff]
                end
            end
        end
    end

    # Set Html_serial as last column
    ret = vcat(out..., cols=:union)
    select!(ret, [filter(x -> x != "Html_serial", names(ret)); "Html_serial"])
    return ret
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
add_MBH_Hb_WuShen2022!(results)
add_MBH_MgII_WuShen2022!(results)
# add_MBH_Ha_ShenLiu2012!(results)
# add_MBH_Ha_Ricci!(results)
add_MBH_Hb_Ricci!(results)
# add_MBH_Ha_L5100_Ricci!(results)
add_MBH_MgII_Ricci!(results)
# add_MBH_Pab_Ricci!(results)
# add_MBH_HeI_Ricci!(results)

results.MBH_mean .= NaN
for i in 1:nrow(results)
    results[i, :MBH_mean] = mean_handle_NaN([results[i, :MBH_Hb_WuShen2022], results[i, :MBH_MgII_WuShen2022]])
end

# Calculates Lbol and Eddington ratios
add_Lbol_eddratio!(results)


# Creates Quality cut columns
results[!, :qcut] = ((results.NPOINTS .> 6000)  .&
                     (results.SNR .> 3))


# Write results in a FITS file
write_fits("$(output_path)/QSFIT_RESULTS.fits", results)
