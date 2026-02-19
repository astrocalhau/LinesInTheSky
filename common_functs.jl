function skip_NaN_missing(v)
    i = findall(.!isnan.(v)  .&  .!ismissing.(v))
    (length(i) == 0)  &&  (return Float64[])
    return v[i]
end



function add_MBH_Hb_WuShen2022!(cc)
    cc[!, :MBH_Hb_WuShen2022] .= NaN
    i = findall((.!isnan.(cc.L5100))          .&
                (cc.QSOcont_reliable .=== 1)  .&
                (cc.Hb_br_reliable   .=== 1))
    cc[i, :MBH_Hb_WuShen2022] .= 0.91 .+ 0.5  .* log10.(cc.L5100[i] .* 1e-2) .+ 2 .* log10.(cc.Hb_br_fwhm[i])
end

function add_MBH_MgII_WuShen2022!(cc)
    cc[!, :MBH_MgII_WuShen2022] .= NaN
    i = findall((.!isnan.(cc.L3000))               .&
                (cc.QSOcont_reliable      .=== 1)  .&
                (cc.MgII_2798_br_reliable .=== 1))
    cc[i, :MBH_MgII_WuShen2022] .= 0.74 .+ 0.62 .* log10.(cc.L3000[i] .* 1e-2) .+ 2 .* log10.(cc.MgII_2798_br_fwhm[i])
end

function add_MBH_Ha_ShenLiu2012!(cc)
    cc[!, :MBH_Ha_ShenLiu2012] .= NaN
    i = findall(cc.Ha_br_reliable .=== 1)
    Ha_norm = cc.Ha_br_norm  # TODO: add narrow component?  also in FWHM?
    cc[i, :MBH_Ha_ShenLiu2012] .= 2.216 .+ 0.564 .* log10.(Ha_norm[i] .* 1e-2) .+ 1.821 .* log10.(cc.Ha_br_fwhm[i])
end


function add_MBH_Ha_Ricci!(cc)
    cc[!, :MBH_Ha_Ricci] .= NaN
    i = findall(cc.Ha_br_reliable .=== 1)
    # Log Mbh/Msun =7.072 +0.563*log(LHalpha/10^44 erg/s)+2log(FWHM/1000 km/s)
    cc[i, :MBH_Ha_Ricci] .= 7.072 .+ 0.563 .* log10.(cc.Ha_br_norm[i] .* 1e-2) .+ 2 .* log10.(cc.Ha_br_fwhm[i] ./ 1000.)
end

function add_MBH_Ha_L5100_Ricci!(cc)
    cc[!, :MBH_Ha_L5100_Ricci] .= NaN
    i = findall((.!isnan.(cc.L5100))          .&
                (cc.QSOcont_reliable .=== 1)  .&
                (cc.Ha_br_reliable .=== 1))
    # 6.779 + 0.650*log10(L5100/1e+44 erg/s) + 2*log10(FWHM_Halpha /1000 km/s)
    cc[i, :MBH_Ha_L5100_Ricci] .= 6.779 .+ 0.65 .* log10.(cc.L5100[i] .* 1e-2) .+ 2 .* log10.(cc.Ha_br_fwhm[i] ./ 1000.)
end


function add_MBH_Hb_Ricci!(cc)
    cc[!, :MBH_Hb_Ricci] .= NaN
    i = findall((.!isnan.(cc.L5100))          .&
                (cc.QSOcont_reliable .=== 1)  .&
                (cc.Hb_br_reliable   .=== 1))
    # Log Mbh/Msun = 6.721 + 0.650*log(L5100/10^44 erg/s) +2log(FWHM/1000 km/s)
    cc[i, :MBH_Hb_Ricci] .= 6.721 .+ 0.65 .* log10.(cc.L5100[i] .* 1e-2) .+ 2 .* log10.(cc.Hb_br_fwhm[i] ./ 1000.)
end

function add_MBH_MgII_Ricci!(cc)
    cc[!, :MBH_MgII_Ricci] .= NaN
    i = findall((.!isnan.(cc.L3000))               .&
                (cc.QSOcont_reliable      .=== 1)  .&
                (cc.MgII_2798_br_reliable .=== 1))
    # Log Mbh/Msun=6.906 + 0.609*log(L3000/10^44 erg/s) +2log(FWHM/1000 km/s)
    cc[i, :MBH_MgII_Ricci] .= 6.906 .+ 0.609 .* log10.(cc.L3000[i] .* 1e-2) .+ 2 .* log10.(cc.MgII_2798_br_fwhm[i] ./ 1000.)
end


function add_MBH_Pab_Ricci!(cc)
    cc[!, :MBH_Pab_Ricci] .= NaN
    i = findall(cc.Pab_br_reliable .=== 1)
    # log(M_BH/M_sun) = 7.94+0.872(+/- 0.040) x { 2x[log(FWHM/1e+4 km/s)] + 0.5x[log(L_Pa_b/1e+40 erg/s]}
    cc[i, :MBH_Pab_Ricci] .= 7.94 .+ 0.872 .* (2 .* log10.(cc.Pab_br_fwhm[i] ./ 1e4) .+ 0.5 .* log10.(cc.Pab_br_norm[i] .* 1e2))
end

function add_MBH_HeI_Ricci!(cc)
    cc[!, :MBH_HeI_Ricci] .= NaN
    i = findall(cc.HeI_10832_br_reliable .=== 1)
    # Log(M_BH/M_sun) = 7.86 + 2x[log(FWHM/1e+4 km/s)] + 0.5x[log(L_HeI erg/s) – 39.55]
    cc[i, :MBH_HeI_Ricci] .= 7.86 .+ 2 .* log10.(cc.HeI_10832_br_fwhm[i] ./ 1e4) .+ 0.5 .* (log10.(cc.HeI_10832_br_norm[i]) .+ 42 .- 39.55)
end

function add_Lbol_eddratio_Euclid!(cc)
    j = findall(skip_NaN_missing((cc.Pab_br_reliable .===1) .&
    				    (cc.Pab_br_norm .> 0 )))
    k = findall(skip_NaN_missing((cc.HeI_10832_br_reliable .===1) .&
    				    (cc.HeI_10832_br_norm .> 0 )))
    L_Pab = log10.(cc.Pab_br_norm[j] .* 1.0E42)
    L_HeI = log10.(cc.HeI_10832_br_norm[k] .* 1.0E42)
    L_sun = 3.828E33 # erg/s
    L_2_10keV_Pab = (10 .^ (L_Pab .+ 2.13)) ./ L_sun #From L_Pab to LX
    L_2_10keV_HeI = (10 .^ (L_HeI .+ 2.13)) ./ L_sun #From L_HeI to LX
    KX_Pab = 15.33 .* (1 .+ (log10.(L_2_10keV_Pab) ./ 11.48) .^ 16.2) # bolometric correction
    KX_HeI = 15.33 .* (1 .+ (log10.(L_2_10keV_HeI) ./ 11.48) .^ 16.2) # bolometric correction
    
    cc.Lbol_3000 = 5.15e44 .* cc.L3000 .* 1e-2
    cc.Lbol_5100 = 9.26e44 .* cc.L5100 .* 1e-2
    cc[!, :Lbol_LX_Pab] .= NaN
    cc[!, :Lbol_LX_HeI] .= NaN
    allowmissing!(cc, :Lbol_LX_Pab)
    allowmissing!(cc, :Lbol_LX_HeI)
    cc[j, :Lbol_LX_Pab] .= (L_2_10keV_Pab .* L_sun) .* KX_Pab
    cc[k, :Lbol_LX_HeI] .= (L_2_10keV_HeI .* L_sun) .* KX_HeI
    cc.Lbol_mean .= 0.
    allowmissing!(cc, :Lbol_mean)
    for i in 1:nrow(cc)
        try
	        cc[i, :Lbol_mean] = mean(skip_NaN_missing([cc[i, :Lbol_3000], cc[i, :Lbol_5100], cc[i, :Lbol_LX_Pab], cc[i, :Lbol_LX_HeI]]))
        catch
            cc[i, :Lbol_mean] = missing
        end
    end

    # Calculates Eddington ratios
    cc.Ledd_mean = 1.26e38 * 10 .^cc.MBH_mean
    cc.Edd_ratio = cc.Lbol_mean ./ cc.Ledd_mean
end

function add_Lbol_eddratio_DESI!(cc)
    cc.Lbol_3000 = 5.15e44 .* cc.L3000 .* 1e-2
    cc.Lbol_5100 = 9.26e44 .* cc.L5100 .* 1e-2
    cc.Lbol_mean .= 0.
    allowmissing!(cc, :Lbol_mean)
    for i in 1:nrow(cc)
        try
	        cc[i, :Lbol_mean] = mean(skip_NaN_missing([cc[i, :Lbol_3000], cc[i, :Lbol_5100]]))
        catch
            cc[i, :Lbol_mean] = missing
        end
    end

    # Calculates Eddington ratios
    cc.Ledd_mean = 1.26e38 * 10 .^cc.MBH_mean
    cc.Edd_ratio = cc.Lbol_mean ./ cc.Ledd_mean
end




