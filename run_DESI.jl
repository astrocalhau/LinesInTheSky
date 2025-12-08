using Revise
using Base.Threads, FITSIO, DataFrames, DataStructures, Printf, Statistics, StatsBase
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer
using SkyCoords, DustExtinction
using LinesInTheSky

include("utils.jl")


function analyze_spec(input_path, output_path, row; clob=false)
    mkpath("$(output_path)/JSON")

    input_filename = "$(input_path)/TXT/$(row[:id_DESI_DR1]).txt"
    isfile(input_filename)  ||  error("Input file $input_filename do not exists.")

    output_filename = "$(output_path)/JSON/$(row[:id_DESI_DR1]).json.gz"
    if isfile(output_filename)  &&  !clob
        println("Output file $output_filename already exists. Skipping.")
        return
    end

    println(); println()
    @info "Analyzing file $(input_filename)"

    gal_coords = convert(GalCoords, ICRSCoords(row[:ra] * pi/180., row[:dec] * pi/180))
    ebv = SFD98Map()(gal_coords.l, gal_coords.b)

    spec = Spectrum(Val(:ASCII), input_filename, columns=[1,2,3], resolution=2857, label=string(row[:id_DESI_DR1]))
    recipe = CRecipe{Type1}(redshift=row[:Z], use_host_template=true, Av=ebv * 3.1)
    res = analyze(recipe, spec)

    QSFit.serialize(output_filename, res, compress=true)
    return res
end

function read_results(output_path, row)
    filename = "$(output_path)/JSON/$(row.id_DESI_DR1).json.gz"
    isfile(filename)  ||  error("File $filename do not exists.")

    @info "Reading $filename ..."
    res = QSFit.deserialize(filename)
    out = OrderedDict(:ID_DESI => row.id_DESI_DR1, :ID_EUCLID => row.object_id,
                      :Redshift => row.Z, :Source => row.CAT, :redchisq => res.fsumm.fitstat,
                      :NPOINTS => res.fsumm.ndata, :SNR => res.post[:Data_stats][:SNR], :DER_SNR => res.post[:Data_stats][:DER_SNR], :nneg => res.post[:Data_stats][:nneg],
                      :L3000 => res.post[:Continuum_luminosity][:l3000],
                      :L5100 => res.post[:Continuum_luminosity][:l5100])
    for (cname, comp) in res.bestfit
        (cname in [:QSOcont, :Galaxy, :Ironuv, :Ironoptbr, :Ironoptna,
                   :Ha_br, :Ha_na, :Hb_br, :Hb_na, :Pab_br, :HeI_10832_br,
                   :MgII_2798_br, :OIII_4959, :OIII_5007, :OIII_5007_bw])  ||  continue

        out[Symbol(cname, :_reliable)] = Int.(!(cname in keys(res.post[:Issues])  &&
            ((:norm in keys(res.post[:Issues][cname]))  ||
                                                (:fwhm in keys(res.post[:Issues][cname])))))
        for (pname, par) in GModelFit.getparams(comp)
            colname = Symbol(cname, :_, pname)
            if isnothing(par.patch)
                out[       colname        ] = par.val
                out[Symbol(colname, :_unc)] = par.unc
            else
                out[      colname        ]  = par.actual
                out[Symbol(colname, :_unc)] = NaN
            end
        end
    end

    for assoc in [:Ha_br_assoc, :Hb_br_assoc]
        if assoc in keys(res.post)
            out[Symbol(assoc, :_norm)] = res.post[assoc][:norm]
            out[Symbol(assoc, :_fwhm)] = res.post[assoc][:fwhm]
            out[Symbol(assoc, :_voff)] = res.post[assoc][:voff]
        end
    end

    return DataFrame(out)
end

include("common_functs.jl")
function calculate_additional_columns!(results)
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
        results[i, :MBH_mean] = mean(skip_NaN_missing([results[i, :MBH_Hb_WuShen2022], results[i, :MBH_MgII_WuShen2022]]))
    end

    # Calculates Lbol and Eddington ratios
    add_Lbol_eddratio!(results)

    # Creates Quality cut columns
    results[!, :good] = Int.((results.NPOINTS .> 6000)                   .&
                             (results.nneg ./ results.NPOINTS .< 0.1)    .&
                             (results.DER_SNR .> 3))
end


function run_DESI()
    input_path  = "input/input_DESI"
    output_path = "results_DESI"

    # Read input catalog
    f = FITS("$(input_path)/Q1_DESI_DR1_QSOAGNVAC_QSO.fits")
    catalog = DataFrame(f[2])
    close(f)

    # Run analysis
    run_analysis(input_path, output_path, catalog)

    # Read results from JSON files
    results = read_all_results(output_path, catalog)

    # Calculate additional columns
    calculate_additional_columns!(results)

    # Write results in a FITS file
    write_fits("$(output_path)/QSFIT_RESULTS.fits", results)

    return input_path, output_path, catalog, results
end

# input_path, output_path, catalog, results = run_DESI()
