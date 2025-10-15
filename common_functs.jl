function mean_handle_NaN(v)
    i = findall(.!isnan.(v)  .&  .!ismissing.(v))
    (length(i) == 0)  &&  (return NaN)
    return mean(v[i])
end

function cont_lambdaLlambda(model, wavelength)
    try
        return Dierckx.Spline1D(coords(domain(model)), model(:QSOcont), k=1, bc="error")(wavelength) * wavelength * 1e-2
    catch
        return NaN
    end
end


function add_MBH_Hb_WuShen2022!(cc)
    cc[!, :MBH_Hb_WuShen2022] .= NaN
    i = findall((.!isnan.(cc.L5100))          .&
                (cc.QSOcont_reliable .=== 1)  .&
                (cc.Hb_br_reliable   .=== 1))
    cc[i, :MBH_Hb_WuShen2022] .= 0.91 .+ 0.5  .* log10.(cc.L5100[i]) .+ 2 .* log10.(cc.Hb_br_fwhm[i])
end

function add_MBH_MgII_WuShen2022!(cc)
    cc[!, :MBH_MgII_WuShen2022] .= NaN
    i = findall((.!isnan.(cc.L3000))               .&
                (cc.QSOcont_reliable      .=== 1)  .&
                (cc.MgII_2798_br_reliable .=== 1))
    cc[i, :MBH_MgII_WuShen2022] .= 0.74 .+ 0.62 .* log10.(cc.L3000[i]) .+ 2 .* log10.(cc.MgII_2798_br_fwhm[i])
end

function add_MBH_Ha_ShenLiu2012!(cc)
    cc[!, :MBH_Ha_ShenLiu2012] .= NaN
    i = findall(cc.Ha_br_reliable .=== 1)
    Ha_norm = cc.Ha_br_norm  # TODO: add narrow component?  also in FWHM?
    cc[i, :MBH_Ha_ShenLiu2012] .= 2.216 .+ 0.564 .* log10.((Ha_norm[i]) .* 1e-2) .+ 1.821 .* log10.(cc.Ha_br_fwhm[i])
end


function add_Lbol_eddratio!(cc)
    cc.Lbol_3000 = 5.15e44 .* cc.L3000
    cc.Lbol_5100 = 9.26e44 .* cc.L5100
    cc.Lbol_mean .= 0.
    allowmissing!(cc, :Lbol_mean)
    for i in 1:nrow(cc)
        try
	        cc[i, :Lbol_mean] = mean_handle_NaN([cc[i, :Lbol_3000], cc[i, :Lbol_5100]])
        catch
            cc[i, :Lbol_mean] = missing
        end
    end

    # Calculates Eddington ratios
    cc.Ledd_mean = 1.26e38 * 10 .^cc.MBH_mean
    cc.Edd_ratio = cc.Lbol_mean ./ cc.Ledd_mean
end


function write_fits(filename, results)
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
    f = FITS(filename, "w")
    write(f, data)
    close(f)
end
