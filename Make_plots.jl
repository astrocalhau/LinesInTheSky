using Statistics, StatsBase, Plots, StatsPlots, LaTeXStrings, DataFrames, FITSIO, Printf, QSFit, JSON, SortMerge

######################################################################
# Script to create the various figures for the Euclid Q1 Paper "Lines in the Sky"

mkpath("Figures")

f = FITS("results_Euclid/QSFIT_RESULTS.fits")
euclid = DataFrame(f[2])
close(f)

f = FITS("results_DESI/QSFIT_RESULTS.fits")
desi = DataFrame(f[2])
close(f)

# Matching DESI and Euclid catalogs
jj = sortmerge(euclid.ID, desi.ID_EUCLID)
# @assert all(countmatch(jj, 2) .== 1)
# TODO: why do we have DESI spectra with no Euclid counterpart????
euclid_match = euclid[jj[1], :]
desi = desi[jj[2], :]

# Identify subsets in Euclid catalog
sub = (FU      = (euclid.FU .== 1),
       DESI    = (euclid.DESI .== 1),
       QUBRICS = (euclid.QUBRICS .== 1),
       lowz    = (        euclid.Redshift .< 0.8),
       cosmo   = (0.8 .<= euclid.Redshift .< 1.9),
       highz   = (1.9 .<= euclid.Redshift),
       hostgal = (euclid.Galaxy_norm .>= 0)) # sources where the host galaxy template was fit (regardless of reliability)

# Identify quality cuts in both Euclid and DESI catalog
function quality_cuts(df)
    out = (good = ((df.good .== 1)  .&  (df.QSOcont_reliable .== 1)),                                       # sources passing the quality cut
           Ha   = ((df.good .== 1)  .&  (df.QSOcont_reliable .== 1)  .&  (df.Ha_br_reliable        .== 1)), # quality cut for Ha-based BH masses
           Hb   = ((df.good .== 1)  .&  (df.QSOcont_reliable .== 1)  .&  (df.Hb_br_reliable        .== 1)), # quality cut for Hb-based BH masses
           MgII = ((df.good .== 1)  .&  (df.QSOcont_reliable .== 1)  .&  (df.MgII_2798_br_reliable .== 1))) # quality cut for MgII-based BH masses
    if "Pab_br_reliable" in names(df)
        out = (out...,
           Pab  = ((df.good .== 1)  .&  (df.QSOcont_reliable .== 1)  .&  (df.Pab_br_reliable       .== 1))) # quality cut for Pab-based BH masses
    end
    if "HeI_10832_br_reliable" in names(df)
        out = (out...,
           HeI  = ((df.good .== 1)  .&  (df.QSOcont_reliable .== 1)  .&  (df.HeI_10832_br_reliable .== 1))) # quality cut for HeI-based BH masses
    end
    return out
end
    
qc = quality_cuts(euclid)              # Quality cuts for the whole Euclid catalog
qc_match = quality_cuts(euclid_match)  # Quality cuts for the subset in Euclid catalog having a DESI spectrum
qc_desi = quality_cuts(desi)           # Quality cuts for the DESI catalog

# General plot settings
GEN_OPTS = (fontfamily="Computer Modern", framestyle=:box, grid=false,
            background_color_legend=nothing, foreground_color_legend=nothing, # , legend_column=-1
            palette=:darkrainbow, # :seaborn_dark6, :rainbow 
            thickness_scaling=1.2)
# , titlefontsize=8, guidefontsize=8, tickfontsize=8, legendfontsize=7
# palette = [:darkred, :darkgreen, :darkblue, :gray, :purple]

HISTO_OPTS = (alpha=1, fillalpha=0.7, fill=true)
 
xfraction(f) = xlims()[1] + (xlims()[2] - xlims()[1]) * f
yfraction(f) = ylims()[1] + (ylims()[2] - ylims()[1]) * f

function add_μσ!(vv; x=0.03, y1=0.8, y2=0.6)
    μ = median(vv)
    σ = mad(vv, normalize=true)
    sμ = @sprintf("%.2f", μ)
    sσ = @sprintf("%.2f", σ)
    vline!([median(vv)]              , label="", color=:black, linewidth=3)
    vline!( median(vv) .+ [1,-1] .* σ, label="", color=:black, linewidth=2, linestyle=:dash)
    annotate!([xfraction(x)], [yfraction(y1)], text(L"$\tilde{\mu}=%$sμ $"   , 10, :black, :left))
    annotate!([xfraction(x)], [yfraction(y2)], text(L"$\tilde{\sigma}=%$sσ $", 10, :black, :left))
end


######################################################################
# Redshift histograms
let
    plot(; GEN_OPTS..., title=L"$\mathrm{\textbf{Full\ sample}}$", xlabel=L"$z$", ylabel=L"$\mathrm{Counts}$")
    SERIES_OPTS = (bins=0:0.25:maximum(euclid.Redshift), HISTO_OPTS...)
    stephist!(euclid[:          , :Redshift], label=L"$\mathrm{Merged\ sample}$"          ; SERIES_OPTS...)
    stephist!(euclid[sub.FU     , :Redshift], label=L"$\mathrm{Fu\ et\ al.\ (in\ prep.)}$"; SERIES_OPTS...)
    stephist!(euclid[sub.DESI   , :Redshift], label=L"$\mathrm{DESI}$"                    ; SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS, :Redshift], label=L"$\mathrm{QUBRICS}$"                 ; SERIES_OPTS...)
    savefig("Figures/Redshift_Hist_Figure.pdf")

    plot(; GEN_OPTS..., title=L"$\mathrm{\textbf{Good\ quality\ sample}}$", xlabel=L"$z$", ylabel=L"$\mathrm{Counts}$")
    stephist!(euclid[                 qc.good, :Redshift], label=L"$\mathrm{Merged\ sample}$"          ; SERIES_OPTS...)
    stephist!(euclid[sub.FU       .&  qc.good, :Redshift], label=L"$\mathrm{Fu\ et\ al.\ (in\ prep.)}$"; SERIES_OPTS...)
    stephist!(euclid[sub.DESI     .&  qc.good, :Redshift], label=L"$\mathrm{DESI}$"                    ; SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS  .&  qc.good, :Redshift], label=L"$\mathrm{QUBRICS}$"                 ; SERIES_OPTS...)
    savefig("Figures/Redshift_Hist_Quality_Figure.pdf")
end



######################################################################
# BH mass histograms
let
    SERIES_OPTS = (bins=6.:0.5:10.5, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm)) # <-- these are necessary to reduce space between the subplots
    accum = Vector{Any}()
    vv = filter(!isnan, euclid[qc.MgII, :MBH_MgII_WuShen2022]); push!(accum, stephist(vv, label=L"$\mathrm{MgII}$"   , color=length(accum)+1; SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.Hb  , :MBH_Hb_WuShen2022])  ; push!(accum, stephist(vv, label=L"$\mathrm{H}\beta$" , color=length(accum)+1; SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.Ha  , :MBH_Ha_ShenLiu2012]) ; push!(accum, stephist(vv, label=L"$\mathrm{H}\alpha$", color=length(accum)+1; SERIES_OPTS..., xformatter=_->"", ylabel=L"$\mathrm{Counts}$")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.HeI , :MBH_HeI_Ricci])      ; push!(accum, stephist(vv, label=L"$\mathrm{He\,I}$"  , color=length(accum)+1; SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.Pab , :MBH_Pab_Ricci])      ; push!(accum, stephist(vv, label=L"$\mathrm{Pa}\beta$", color=length(accum)+1; SERIES_OPTS..., xlabel=L"$\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})$", bottom_margin=(-1., :mm))); add_μσ!(vv)
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS..., legend=:topright)
    savefig("Figures/BH_mass.pdf")

    # Comparison between Ha and Hb
    plot(; GEN_OPTS..., xlabel=L"$\log_{10}(M_{\mathrm{H}\alpha}/M_{\mathrm{H}\beta})$", ylabel=L"$\mathrm{Counts}$")
    ii = findall(qc.Ha  .&  qc.Hb)
    vv = filter(!isnan, euclid[ii, :MBH_Ha_ShenLiu2012] .- euclid[ii, :MBH_Hb_WuShen2022])
    stephist!(vv, label=L"$\mathrm{H}\alpha\ vs\ \mathrm{H}\beta$", bins=minimum(vv):0.25:maximum(vv); HISTO_OPTS...)
    add_μσ!(vv, y1=0.8, y2=0.7)
    savefig("Figures/BH_mass_Ha_vs_Hb.pdf")
end


######################################################################
# QSO Continuum alpha index histogram
let
    function add_series!(vv, label, color, SERIES_OPTS; showmean=true, kws...)
        μ = median(vv)
        σ = mad(vv, normalize=true)
        sμ = @sprintf("%.2f", μ)
        sσ = @sprintf("%.2f", σ)
        if showmean
            stephist!(vv, label=label * L"\ (\tilde{\mu}=%$sμ)", color=color; SERIES_OPTS..., kws...)
            vline!([median(vv)], label="", color=color, linewidth=2, linestyle=:dash)
        else
            stephist!(vv, label=label, color=color; SERIES_OPTS..., kws...)
        end
    end
    
    plot(; GEN_OPTS..., legend=:topright, xlabel=L"$\alpha_{\lambda}$", ylabel=L"$\mathrm{Counts}$")
    SERIES_OPTS = (bins=minimum(euclid.QSOcont_alpha):0.25:maximum(euclid.QSOcont_alpha), HISTO_OPTS...)
    add_series!(euclid[qc.good              , :QSOcont_alpha], L"$\mathrm{Merged\ sample}$", :black, SERIES_OPTS, showmean=false, fill=false, linewidth=4)
    add_series!(euclid[qc.good .&  sub.lowz , :QSOcont_alpha], L"$z < 0.8$"                , 1     , SERIES_OPTS)
    add_series!(euclid[qc.good .&  sub.cosmo, :QSOcont_alpha], L"$0.8 < z < 1.9$"          , 2     , SERIES_OPTS)
    add_series!(euclid[qc.good .&  sub.highz, :QSOcont_alpha], L"$z > 1.9$"                , 3     , SERIES_OPTS)
    add_series!(desi[qc_desi.good           , :QSOcont_alpha], L"DESI"                     , 4     , SERIES_OPTS, z_order=2)
    savefig("Figures/QSOcont_alpha_Hist_Redshift_2.pdf")
end


######################################################################
# Chi2 vs SNR plots
let
    plot(; GEN_OPTS..., xlabel=L"$\mathrm{Reduced}\ \chi^2$", ylabel=L"$\mathrm{Spectrum\ S/N}$")
    scatter!(euclid[:      , :redchisq], euclid[:      , :DER_SNR], label=L"$\mathrm{Full\ sample}$", ms=2, markerstrokewidth=0)
    scatter!(euclid[qc.good, :redchisq], euclid[qc.good, :DER_SNR], label=L"$\mathrm{Good\ sample}$", ms=4, ma=0.6, markershape=:circ)
    plot!(xaxis=:log10, yaxis=:log10)
    savefig("Figures/Chi2vsSNR_cut_Fig.pdf")
end


######################################################################
# H magnitude histogram
let
    SERIES_OPTS = (bins=minimum(euclid.Hmag):0.25:maximum(euclid.Hmag), HISTO_OPTS...)
    plot(; GEN_OPTS..., xlabel=L"$H_E\ \mathrm{magnitude}$", ylabel=L"$\mathrm{Counts}$", title=L"$\mathrm{\textbf{Good\ sample}}$", legend=:topleft)
    stephist!(euclid[               qc.good, :Hmag], label=L"$\mathrm{Merged\ sample}$"          ; SERIES_OPTS...)
    stephist!(euclid[sub.FU      .& qc.good, :Hmag], label=L"$\mathrm{Fu\ et\ al.\ (in\ prep.)}$"; SERIES_OPTS...)
    stephist!(euclid[sub.DESI    .& qc.good, :Hmag], label=L"$\mathrm{DESI}$"                    ; SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS .& qc.good, :Hmag], label=L"$\mathrm{QUBRICS}$"                 ; SERIES_OPTS...)
    savefig("Figures/Hmag_Hist_source_Fig.pdf")
end


######################################################################
# Hmag vs redshift
let
    SERIES_OPTS = (ms=6, ma=0.6)
    plot(; GEN_OPTS..., title=L"$\mathrm{\textbf{Good\ sample}}$", xlabel=L"$z$", ylabel=L"$H_E$", ylims=(14, 23))
    i = findall(sub.FU      .& qc.good); scatter!(euclid[sub.FU      .& qc.good, :Redshift], euclid[i, :Hmag], label=L"$\mathrm{Fu\ et\ al.}$", markershape=:rect     ; SERIES_OPTS...)
    i = findall(sub.DESI    .& qc.good); scatter!(euclid[sub.DESI    .& qc.good, :Redshift], euclid[i, :Hmag], label=L"$\mathrm{DESI}$"       , markershape=:pentagon ; SERIES_OPTS...)
    i = findall(sub.QUBRICS .& qc.good); scatter!(euclid[sub.QUBRICS .& qc.good, :Redshift], euclid[i, :Hmag], label=L"$\mathrm{QUBRICS}$"    , markershape=:utriangle; SERIES_OPTS...)
    savefig("Figures/Hmag_vs_z_cut_Fig.pdf")
end


######################################################################
# Emission line histograms
let
    SERIES_OPTS = (bins=40.5:0.25:45, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm)) # <-- these are necessary to reduce space between the subplots
    accum = Vector{Any}()
    push!(accum, stephist(log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42, label=L"$\mathrm{Mg\,II}$";  color=length(accum)+1, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Hb  , :Hb_br_norm])        .+ 42, label=L"$\mathrm{H\beta}$";  color=length(accum)+1, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Ha  , :Ha_br_norm])        .+ 42, label=L"$\mathrm{H\alpha}$"; color=length(accum)+1, SERIES_OPTS..., xformatter=_->"", ylabel=L"$\mathrm{Counts}$"))
    push!(accum, stephist(log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42, label=L"$\mathrm{He\,I}$";   color=length(accum)+1, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Pab , :Pab_br_norm])       .+ 42, label=L"$\mathrm{Pa}\beta$"; color=length(accum)+1, SERIES_OPTS..., xlabel=L"$\log_{10}(L_{\mathrm{line}}/\mathrm{erg\ s^{-1}})$", bottom_margin=(-1., :mm)))
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS...)
    savefig("Figures/Lines_lum.pdf")
end

let
    SERIES_OPTS = (bins=2.5:0.2:5, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm)) # <-- these are necessary to reduce space between the subplots
    accum = Vector{Any}()
    push!(accum, stephist(log10.(euclid[qc.MgII, :MgII_2798_br_fwhm]), label=L"$\mathrm{Mg\,II}$";  color=1, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Hb  , :Hb_br_fwhm])       , label=L"$\mathrm{H\beta}$";  color=2, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Ha  , :Ha_br_fwhm])       , label=L"$\mathrm{H\alpha}$"; color=3, SERIES_OPTS..., xformatter=_->"", ylabel=L"$\mathrm{Counts}$"))
    push!(accum, stephist(log10.(euclid[qc.HeI , :HeI_10832_br_fwhm]), label=L"$\mathrm{He\,I}$";   color=4, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Pab , :Pab_br_fwhm])      , label=L"$\mathrm{Pa}\beta$"; color=5, SERIES_OPTS..., xlabel=L"$\log_{10}(\mathrm{FWHM/km\ s^{-1}})$", bottom_margin=(-1., :mm)))

    SERIES_OPTS = (SERIES_OPTS..., bins=-500:200.:500)
    push!(accum, stephist(       euclid[qc.MgII, :MgII_2798_br_voff] , label=L"$\mathrm{Mg\,II}$";  color=1, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.Hb  , :Hb_br_voff]        , label=L"$\mathrm{H\beta}$";  color=2, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.Ha  , :Ha_br_voff]        , label=L"$\mathrm{H\alpha}$"; color=3, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.HeI , :HeI_10832_br_voff] , label=L"$\mathrm{He\,I}$";   color=4, SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.Pab , :Pab_br_voff]       , label=L"$\mathrm{Pa}\beta$"; color=5, SERIES_OPTS..., xlabel=L"$\mathrm{V_{off}/km\ s^{-1}}$", bottom_margin=(-1., :mm)))

    accum = permutedims(reshape(accum,     div(length(accum), 2), 2)) # reorder subplots so that they appear in the correct order
    plot(reshape(accum, :)..., layout=grid(div(length(accum), 2), 2); GEN_OPTS...)
    savefig("Figures/Lines_FWHM_Voff.pdf")
end


######################################################################
# Bol Lum Mean Histogram
let
    plot(; GEN_OPTS..., legend=:topleft, xlabel=L"$\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})$", ylabel=L"$\mathrm{Counts}$")
    vv = filter(!isnan, euclid[qc.good             , :MBH_mean]);  stephist!(vv, label=L"$\mathrm{Good\ sample}$"; HISTO_OPTS...)
    vv = filter(!isnan, euclid[qc.good .& sub.cosmo, :MBH_mean]);  stephist!(vv, label=L"$0.8 < z < 1.9$"        ; HISTO_OPTS...)
    savefig("Figures/Mean_BH_Mass.pdf")
end

let
    plot(; GEN_OPTS..., legend=:topleft, xlabel=L"$\log_{10}(L_{bol}/\mathrm{erg\ s^{-1}})$", ylabel=L"$\mathrm{Counts}$")
    vv = filter(!isnan, euclid[qc.good             , :Lbol_mean]); stephist!(log10.(vv), label=L"$\mathrm{Good\ sample}$"; HISTO_OPTS...)
    vv = filter(!isnan, euclid[qc.good .& sub.cosmo, :Lbol_mean]); stephist!(log10.(vv), label=L"$0.8 < z < 1.9$"        ; HISTO_OPTS...)
    savefig("Figures/Mean_Lbol.pdf")
end


######################################################################
# DESI BH Mass - Euclid Ha BH mass
let
    plot(; GEN_OPTS..., xlabel=L"$\log_{10}(M_{\mathrm{BH,\ Euclid,\ H\alpha}} / M_{\mathrm{BH,\ DESI,\ MgII}})$", ylabel=L"$\mathrm{Counts}$")
    ii = qc_match.Ha .& qc_desi.MgII
    vv = filter(!isnan, euclid_match[ii, :MBH_Ha_ShenLiu2012] .- desi[ii, :MBH_MgII_WuShen2022])
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("Figures/BH_mass_cmpDESI_histo.pdf")

    plot(; GEN_OPTS..., xlabel=L"$\log_{10}(M_{\mathrm{BH,\ Euclid,\ H\alpha}}$", ylabel=L"M_{\mathrm{BH,\ DESI,\ MgII}})", xlims=(7.5, 10), ylims=(7.5,10))
    ii = qc_match.Ha .& qc_desi.MgII
    scatter!(euclid_match[ii, :MBH_Ha_ShenLiu2012], desi[ii, :MBH_MgII_WuShen2022], label="")
    plot!([xlims()...], [ylims()...], label = L"$1:1$", linecolor=:black, ls=:dash, lw=2)
    savefig("Figures/BH_mass_cmpDESI_scatter.pdf")
end

let
    plot(; GEN_OPTS..., xlabel=L"$\log_{10}(\mathrm{L}_{\mathrm{Euclid,\ H\alpha}}\ /\ \mathrm{L}_{\mathrm{DESI,\ MgII}})$", ylabel=L"$\mathrm{Counts}$")
    ii = qc_match.Ha .& qc_desi.MgII
    vv = log10.(filter(!isnan, euclid_match[ii, :Ha_br_norm] ./ desi[ii, :MgII_2798_br_norm]))
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("Figures/Line_ratio_cmpDESI.pdf")
end


######################################################################
# Chi2 vs SNR plots
let
    plot(; GEN_OPTS..., xlabel=L"$\mathrm{Reduced}\ \chi^2$", ylabel=L"$\mathrm{Spectrum\ S/N}$")
    scatter!(desi[:           , :redchisq], desi[:           , :DER_SNR], label=L"$\mathrm{Full\ sample}$", ms=2, markerstrokewidth=0)
    scatter!(desi[qc_desi.good, :redchisq], desi[qc_desi.good, :DER_SNR], label=L"$\mathrm{Good\ sample}$", ms=4, ma=0.6, markershape=:circ)
    plot!(xaxis=:log10, yaxis=:log10)
    savefig("Figures/Chi2vsSNR_cut_Fig_DESI.pdf")
end


#############################################################################
#			SDSS plots 				#

Data = FITS("results_Euclid/QSFIT_RESULTS.fits")
df = DataFrame(Data[2])
close(Data)

#Quality cut
qualcutHa = df[(df.good .== 1) .&& (df.Ha_br_reliable .==1) .&& coalesce.(df.Ha_br_fwhm.>2000, false) .&& coalesce.(df.Ha_br_fwhm.<10000, false), :]
qualcutHb = df[(df.good .== 1) .&& (df.Hb_br_reliable .==1) .&& coalesce.(df.Hb_br_fwhm.>2000, false) .&& coalesce.(df.Hb_br_fwhm.<10000, false), :]
qualcutMgII = df[(df.good .== 1) .&& (df.MgII_2798_br_reliable .==1) .&& coalesce.(df.MgII_2798_br_fwhm.>2000 .&& coalesce.(df.MgII_2798_br_fwhm.<10000, false), false), :] #coalesce needed for handling "missing" values
qualcutPab = df[(df.good .== 1) .&& (df.Pab_br_reliable .==1) .&& coalesce.(df.Pab_br_fwhm.>2000 .&& coalesce.(df.Pab_br_fwhm.<10000, false), false), :] #coalesce needed for handling "missing" values
qualcutHeI = df[(df.good .== 1) .&& (df.HeI_10832_br_reliable .==1) .&& coalesce.(df.HeI_10832_br_fwhm.>2000 .&& coalesce.(df.HeI_10832_br_fwhm.<10000, false), false), :] #coalesce needed for handling "missing" values

########################################################################################################
########################################################################################################

#########################################################################################################
#########################################################################################################

MBHPabZcut = @df qualcutPab scatter(:Redshift, :MBH_Pab_Ricci, label=L"$\mathrm{Pa\beta}$", mc=:azure1, markershape=:dtriangle, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, legend=:bottomright,framestyle=:box)

MBHHeIZcut = @df qualcutHeI scatter!(:Redshift, :MBH_HeI_Ricci, label=L"$\mathrm{He\,I}$", mc=:skyblue2, markershape=:circle, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, legend=:bottomright,framestyle=:box)


MBHHaZcut = @df qualcutHa scatter!(:Redshift, :MBH_Ha_ShenLiu2012, mc=:skyblue3, markershape=:rect, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, label=L"$\mathrm{H\alpha}$")


MBHHbZcut = @df qualcutHb scatter!(:Redshift, :MBH_Hb_WuShen2022, mc=:royalblue3, markershape=:pentagon, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, label=L"$\mathrm{H\beta}$")


MBHMgIIZcut = @df qualcutMgII scatter!(:Redshift, :MBH_MgII_WuShen2022, label=L"$\mathrm{Mg\,II}$", mc=:black, markershape=:utriangle, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, legend=:bottomright,framestyle=:box)

plot!(formatter=:latex, guidefontsize=18, tickfontsize=18, legendfontsize=18)
xlabel!(L"$z$")
ylabel!(L"$\mathrm{\log_{10}(M_{\mathrm{BH}}/{\mathrm{M_{\odot}}})}$")


#########################################################################################################
#########################################################################################################

#For the Eddington ratio lines
eddration1 = 0.1
#eddration2 = 0.3
eddration3 = 1.
#eddration4 = 5.
eddration5 = 10.

x= range(40,50, length=100)

Edlines1 = log10.((10 .^(x))/(eddration1 * (1.26E38)))
#Edlines2 = log10.((10 .^(x))/(eddration2 * (1.26E38)))
Edlines3 = log10.((10 .^(x))/(eddration3 * (1.26E38)))
#Edlines4 = log10.((10 .^(x))/(eddration4 * (1.26E38)))
Edlines5 = log10.((10 .^(x))/(eddration5 * (1.26E38)))



Edrats1 = plot(x, Edlines1, linecolor=:black, label="", ls=:dash, lw=2)
Edrats2 = plot!(x, Edlines3, linecolor=:green, label="", ls=:dash, lw=2)
Edrats3 = plot!(x, Edlines5, linecolor= :orchid, label="", ls=:dash, lw=2)


MBHHaZcutbol = @df qualcutHa scatter!(log10.(:Lbol_mean), :MBH_Ha_ShenLiu2012, mc=:skyblue3, markershape=:rect, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, label=L"$\mathrm{H\alpha}$")

MBHHbZcutbol = @df qualcutHb scatter!(log10.(:Lbol_mean), :MBH_Hb_WuShen2022, mc=:royalblue3, markershape=:pentagon, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, label=L"$\mathrm{H\beta}$")

MBHMgIIZcutbol = @df qualcutMgII scatter!(log10.(:Lbol_mean), :MBH_MgII_WuShen2022, label=L"$\mathrm{Mg\,II}$", mc=:black, markershape=:utriangle, ms=6, ma=1, dpi=300, guidefontsize=16, tickfontsize=16, legendfontsize=16, grid=false, legend=:bottomright,framestyle=:box)

MBHPab = @df qualcutPab scatter!(log10.(:Lbol_mean), :MBH_Pab_Ricci, label=L"$\mathrm{Pa\beta}$", mc=:azure1, markershape=:dtriangle, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, legend=:bottomright,framestyle=:box)

MBHHeI = @df qualcutHeI scatter!(log10.(:Lbol_mean), :MBH_HeI_Ricci, label=L"$\mathrm{He\,I}$", mc=:skyblue2, markershape=:circle, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, legend=:bottomright,framestyle=:box)


annotate!([44.4], [6.85], text(L"$\eta_{Edd}=0.1$", 18, :black , rotation=20))
annotate!([44.4], [5.85], text(L"$\eta_{Edd}=1$", 18, :green , rotation=20))
annotate!([44.4], [4.85], text(L"$\eta_{Edd}=10$", 18, :orchid , rotation=20))

plot!(formatter=:latex,guidefontsize=18, tickfontsize=18, legendfontsize=18)
xlabel!(L"$\log_{10}(L_{\mathrm{bol}}/\mathrm{erg\ s^{-1}})$")
ylabel!(L"$\log_{10}(M_{\mathrm{BH}}/{\mathrm{M_{\odot}}})$")
ylims!(4,12)
xlims!(42.5,49)

#savefig(SDSSMgIIbol, "MBHvLbol_w_SDSS.png")
#savefig(SDSSMgIIbol, "MBHvLbol_w_SDSS.pdf")
AllSDSS = plot(MBHPabZcut, MBHHaZcutbol, layout=grid(1, 2, widths=(4/8, 4/8)), size=(1600,600), margin=5*Plots.mm, left_margin=10*Plots.mm, bottom_margin=13*Plots.mm, titlefontsize=22, guidefontsize=22, tickfontsize=22, legendfontsize=20, annotationfontsize=20)
plot!(formatter=:latex)
savefig(AllSDSS, "Figures/MBH_SDSS_all.png")
