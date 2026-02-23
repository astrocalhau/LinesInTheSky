using Revise
using Base.Threads, FITSIO, DataFrames, DataStructures, Printf, Statistics, StatsBase, Unitful, TypedJSON
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer, SortMerge
using SkyCoords, DustExtinction
using LinesInTheSky

include("utils.jl")

import QSFit.Spectrum
function Spectrum(::Val{:EUCLID}, file::AbstractString; ndrop=10, resolution=450., kws...)
    f = FITS(file)
    scale = read_header(f[2])["FSCALE"]
    wl   = 1. .* float.(read(f[2], "WAVELENGTH"))
    flux = 1. .* float.(read(f[2], "SIGNAL")) .*scale
    var  = 1. .* float.(read(f[2], "VAR")) .* (scale) .*(scale)
    mask = read(f[2], "MASK")
    close(f)

    i = findall(var .> 0)
    wl   = wl[i]
    flux = flux[i]
    var  = var[i]
    mask = mask[i]

    good = convert(Vector{Bool}, (((mask .== 0) .| (mask .== 2))  .&
                                  (var .> 0)    .&
                                  (flux .> 0)))
    if ndrop > 0
        good[1:ndrop] .= false
        good[end-ndrop+1:end] .= false
    end

    out = Spectrum(wl, flux, sqrt.(var);
                   unit_x = u"angstrom",
                   unit_y = u"erg" / u"s" / u"cm"^2 / u"angstrom",
                   good=good, resolution=resolution, kws...)
    return out
end

function analyze_spec(input_path, output_path, row; clob=false)
    mkpath("$(output_path)/JSON")

    input_filename = "$(input_path)/fits/$(row[:object_id]).fits"
    isfile(input_filename)  ||  error("Input file $input_filename do not exists.")

    output_filename = "$(output_path)/JSON/$(row[:object_id]).json.gz"
    if isfile(output_filename)  &&  !clob
        println("Output file $output_filename already exists. Skipping.")
        return
    end

    println(); println()
    @info "Analyzing file $(input_filename)"

    gal_coords = convert(GalCoords, ICRSCoords(row[:ra] * pi/180., row[:dec] * pi/180))
    ebv = SFD98Map()(gal_coords.l, gal_coords.b)

    spec = Spectrum(Val(:EUCLID), input_filename, label=string(row[:object_id]))
    recipe = CRecipe{WP9Type1IR}(redshift=row[:Z], use_host_template=false, Av=ebv * 3.1, n_nuisance=2)
    resNoHost = analyze(recipe, spec)

    recipe.use_host_template = true
    resWithHost = analyze(recipe, spec)


    if !(:QSOcont in keys(resNoHost.post[:Issues]))  &&
       !(:QSOcont in keys(resWithHost.post[:Issues]))
        if resNoHost.fsumm.fitstat < resWithHost.fsumm.fitstat
            res = resNoHost
        else
            res = resWithHost
        end
    else # simply pick the one with no QSOcont issues
        if !(:QSOcont in keys(resNoHost.post[:Issues]))
            res = resNoHost
        else
            res = resWithHost
        end
    end

    TypedJSON.serialize(output_filename, res, compress=true)
    return res
end

function read_results(output_path, row)
    filename = "$(output_path)/JSON/$(row.object_id).json.gz"
    isfile(filename)  ||  error("File filename do not exists.")

    @info "Reading $filename ..."
    res = TypedJSON.deserialize(filename)
    out = OrderedDict(:ID => row.object_id,
                      :Redshift => row.Z, :Hmag => row.HMAG, :Ref_QUBRICS => row.ref_QUBRICS, :QUBRICS => row.QUBRICS, :DESI => row.DESI, :FU => row.FU, :redchisq => res.fsumm.fitstat,
                      :NPOINTS => res.fsumm.ndata, :SNR => res.post[:Data_stats][:SNR], :DER_SNR => res.post[:Data_stats][:DER_SNR], :nneg => res.post[:Data_stats][:nneg],
                      :L3000 => res.post[:Continuum_luminosity][:l3000],
                      :L5100 => res.post[:Continuum_luminosity][:l5100])
    for (cname, comp) in res.bestfit
        (cname in [:QSOcont, :Galaxy, :Ironuv, :Ironoptbr, :Ironoptna,
                   :Ha_br, :Ha_na, :Hb_br, :Hb_na, :Pab_br, :HeI_10832_br,
                   :MgII_2798_br, :OIII_4959, :OIII_5007, :OIII_5007_bw])  ||  continue

        out[Symbol(cname, :_reliable)] = Int.(!(cname in keys(res.post[:Issues])))
        for (pname, par) in GModelFit.getparams(comp)
            colname = Symbol(cname, :_, pname)
            if isnothing(par.patch)
                out[       colname        ] = par.val
                out[Symbol(colname, :_unc)] = par.unc
            else
                out[       colname        ] = par.actual
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
    add_MBH_Ha_ShenLiu2012!(results)
    add_MBH_Ha_Ricci!(results)
    add_MBH_Hb_Ricci!(results)
    add_MBH_Ha_L5100_Ricci!(results)
    add_MBH_MgII_Ricci!(results)
    add_MBH_Pab_Ricci!(results)
    add_MBH_HeI_Ricci!(results)

    results.MBH_mean .= NaN
    for i in 1:nrow(results)
        results[i, :MBH_mean] = mean(skip_NaN_missing([results[i, :MBH_Hb_WuShen2022], results[i, :MBH_MgII_WuShen2022], results[i, :MBH_Ha_ShenLiu2012],results[i, :MBH_Pab_Ricci],results[i, :MBH_HeI_Ricci]]))
    end

    # Calculates Lbol and Eddington ratios
    add_Lbol_eddratio_Euclid!(results)

    # Creates Quality cut columns
    results[!, :good] = Int.((results.NPOINTS .> 450)                                        .&
                             (results.nneg ./ results.NPOINTS .< 0.1)                        .&
                             (results.DER_SNR .> 3)                                          .&
                             (results.redchisq .< 6)                                         .&
                             (results.QSOcont_reliable .== 1))
end


function run_Euclid()
    input_path  = "input_Euclid"
    output_path = "results_Euclid"

    # Read input catalog
    f = FITS("$(input_path)/catalog.fits")
    catalog = DataFrame(f[2])
    close(f)

    # Use redshifts from DESI whenever available
    f = FITS("input_DESI/Q1_DESI_DR1_QSOAGNVAC_QSO.fits")
    desi = DataFrame(f[2])
    close(f)
    jj = sortmerge(catalog.object_id, desi.object_id)
    catalog[jj[1], :Z] .= desi[jj[2], :Z]

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

# input_path, output_path, catalog, results = run_Euclid()
