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
    @assert all(df[df.good .== 1, :QSOcont_reliable] .== 1)  # QSOcont_reliable is a necessary condition for the source to be a "good" one
    out = (good = ((df.good .== 1)),                                       # sources passing the quality cut
           Ha   = ((df.good .== 1)  .&  (df.Ha_br_reliable        .== 1)), # quality cut for Ha-based BH masses
           Hb   = ((df.good .== 1)  .&  (df.Hb_br_reliable        .== 1)), # quality cut for Hb-based BH masses
           MgII = ((df.good .== 1)  .&  (df.MgII_2798_br_reliable .== 1))) # quality cut for MgII-based BH masses
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
            # palette=:darkrainbow, # :seaborn_dark6, :rainbow
            palette = [:darkred, :darkgreen, :darkblue, :darkorange3, :purple],
            thickness_scaling=1.4)
# , titlefontsize=8, guidefontsize=8, tickfontsize=8, legendfontsize=7

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
    annotate!([xfraction(x)], [yfraction(y1)], text(L"\tilde{\mu}=%$sμ "   , 10, :black, :left))
    annotate!([xfraction(x)], [yfraction(y2)], text(L"\tilde{\sigma}=%$sσ ", 10, :black, :left))
end

ith_color(i) = palette(GEN_OPTS.palette)[i]


######################################################################
# Redshift histograms
let
    plot(; GEN_OPTS..., title="", xlabel=L"z", ylabel=L"\mathrm{Counts}")
    SERIES_OPTS = (bins=0:0.25:maximum(euclid.Redshift), HISTO_OPTS...)
    stephist!(euclid[:          , :Redshift], label=L"\mathrm{Full\ sample}"          ; SERIES_OPTS...)
    stephist!(euclid[sub.FU     , :Redshift], label=L"\mathrm{Fu\ et\ al.\ (in\ prep.)}"; SERIES_OPTS...)
    stephist!(euclid[sub.DESI   , :Redshift], label=L"\mathrm{DESI}"                    ; SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS, :Redshift], label=L"\mathrm{QUBRICS}"                 ; SERIES_OPTS...)
    savefig("Figures/Redshift_Hist_full_sample.pdf")

    plot(; GEN_OPTS..., title="", xlabel=L"z", ylabel=L"\mathrm{Counts}")
    stephist!(euclid[                 qc.good, :Redshift], label=L"\mathrm{Good\ sample}"          ; SERIES_OPTS...)
    stephist!(euclid[sub.FU       .&  qc.good, :Redshift], label=L"\mathrm{Fu\ et\ al.\ (in\ prep.)}"; SERIES_OPTS...)
    stephist!(euclid[sub.DESI     .&  qc.good, :Redshift], label=L"\mathrm{DESI}"                    ; SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS  .&  qc.good, :Redshift], label=L"\mathrm{QUBRICS}"                 ; SERIES_OPTS...)
    savefig("Figures/Redshift_Hist_Quality_sample.pdf")
end


######################################################################
# BH mass histograms
let
    SERIES_OPTS = (bins=6.:0.5:10.5, HISTO_OPTS...,
                   bottom_margin=(-2.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-3.5, :mm), # <-- to reduce wasted space to the left of the plot
                   size=(450, 550)) # <-- Controls the size of the overall plot/image (in pixels only)
    accum = Vector{Any}()
    vv = filter(!isnan, euclid[qc.MgII, :MBH_MgII_WuShen2022]); push!(accum, stephist(vv, label=L"\mathrm{MgII}"   , color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.Hb  , :MBH_Hb_WuShen2022])  ; push!(accum, stephist(vv, label=L"\mathrm{H}\beta" , color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.Ha  , :MBH_Ha_ShenLiu2012]) ; push!(accum, stephist(vv, label=L"\mathrm{H}\alpha", color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"", ylabel=L"\mathrm{Counts}")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.HeI , :MBH_HeI_Ricci])      ; push!(accum, stephist(vv, label=L"\mathrm{He\,I}"  , color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv)
    vv = filter(!isnan, euclid[qc.Pab , :MBH_Pab_Ricci])      ; push!(accum, stephist(vv, label=L"\mathrm{Pa}\beta", color=ith_color(length(accum)+1); SERIES_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})", bottom_margin=(-1., :mm))); add_μσ!(vv)
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS..., legend=:topright)
    savefig("Figures/BH_mass.pdf")

    # Comparison between Ha and Hb
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{H}\alpha}/M_{\mathrm{H}\beta})", ylabel=L"\mathrm{Counts}")
    ii = findall(qc.Ha  .&  qc.Hb)
    vv = filter(!isnan, euclid[ii, :MBH_Ha_ShenLiu2012] .- euclid[ii, :MBH_Hb_WuShen2022])
    stephist!(vv, label=L"\mathrm{H}\alpha\ vs\ \mathrm{H}\beta", bins=minimum(vv):0.25:maximum(vv); HISTO_OPTS...)
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
            out = stephist(vv, label=label * L"\ (\tilde{\mu}=%$sμ)", color=color; SERIES_OPTS..., kws...)
            vline!([median(vv)], label="", color=color, linewidth=2, linestyle=:dash)
        else
            out = stephist(vv, label=label, color=color; SERIES_OPTS..., kws...)
        end
    end
    
    # plot(; GEN_OPTS..., legend=:topright, xlabel=L"\alpha_{\lambda}", ylabel=L"\mathrm{Counts}")
    SERIES_OPTS = (bins=minimum(euclid.QSOcont_alpha):0.25:maximum(euclid.QSOcont_alpha), HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-3.5, :mm), # <-- to reduce wasted space to the left of the plot
                   size=(500, 850)) # <-- Controls size of overall plot
    accum = Vector{Any}()
    push!(accum, add_series!(euclid[qc.good              , :QSOcont_alpha], L"\mathrm{Good\ sample}", :black      , SERIES_OPTS, xformatter=_->"", showmean=false, fill=false, linewidth=4))
    push!(accum, add_series!(euclid[qc.good .&  sub.lowz , :QSOcont_alpha], L"z < 0.8"                , ith_color(1), SERIES_OPTS, xformatter=_->""))
    push!(accum, add_series!(euclid[qc.good .&  sub.cosmo, :QSOcont_alpha], L"0.8 < z < 1.9"          , ith_color(2), SERIES_OPTS, xformatter=_->"", ylabel=L"\mathrm{Counts}"))
    push!(accum, add_series!(euclid[qc.good .&  sub.highz, :QSOcont_alpha], L"z > 1.9"                , ith_color(3), SERIES_OPTS, xformatter=_->""))
    push!(accum, add_series!(desi[qc_desi.good           , :QSOcont_alpha], L"DESI"                   , ith_color(4), SERIES_OPTS, xlabel=L"\alpha_{\lambda}", bottom_margin=(-6., :mm)))
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS..., legend=:topright)
    savefig("Figures/QSOcont_alpha_Hist_Redshift.pdf")
end


######################################################################
# Chi2 vs SNR plots
let
    plot(; GEN_OPTS..., xlabel=L"\mathrm{Reduced}\ \chi^2", ylabel=L"\mathrm{Spectrum\ S/N}")
    scatter!(euclid[:      , :redchisq], euclid[:      , :DER_SNR], label=L"\mathrm{Full\ sample}", ms=2, markerstrokewidth=0)
    scatter!(euclid[qc.good, :redchisq], euclid[qc.good, :DER_SNR], label=L"\mathrm{Good\ sample}", ms=4, ma=0.6, markershape=:circ)
    plot!(xaxis=:log10, yaxis=:log10)
    savefig("Figures/Chi2vsSNR_cut.pdf")
end


######################################################################
# H magnitude histogram
let
    SERIES_OPTS = (bins=minimum(euclid.Hmag):0.25:maximum(euclid.Hmag), HISTO_OPTS...)
    plot(; GEN_OPTS..., xlabel=L"H_E\ \mathrm{magnitude}", ylabel=L"\mathrm{Counts}", title="", legend=:topleft)
    stephist!(euclid[               qc.good, :Hmag], label=L"\mathrm{Good\ sample}"          ; SERIES_OPTS...)
    stephist!(euclid[sub.FU      .& qc.good, :Hmag], label=L"\mathrm{Fu\ et\ al.\ (in\ prep.)}"; SERIES_OPTS...)
    stephist!(euclid[sub.DESI    .& qc.good, :Hmag], label=L"\mathrm{DESI}"                    ; SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS .& qc.good, :Hmag], label=L"\mathrm{QUBRICS}"                 ; SERIES_OPTS...)
    savefig("Figures/Hmag_Hist.pdf")
end


######################################################################
# Hmag vs redshift
let
    SERIES_OPTS = (ms=6, ma=0.6, legend=:bottomright)
    plot(; GEN_OPTS..., title=L"\mathrm{\textbf{Good\ sample}}", xlabel=L"z", ylabel=L"H_E", ylims=(14, 23))
    i = findall(sub.FU      .& qc.good); scatter!(euclid[sub.FU      .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{Fu\ et\ al.}", markershape=:rect     ; SERIES_OPTS...)
    i = findall(sub.DESI    .& qc.good); scatter!(euclid[sub.DESI    .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{DESI}"       , markershape=:pentagon ; SERIES_OPTS...)
    i = findall(sub.QUBRICS .& qc.good); scatter!(euclid[sub.QUBRICS .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{QUBRICS}"    , markershape=:utriangle; SERIES_OPTS...)
    savefig("Figures/Hmag_vs_z_cut.pdf")
end


######################################################################
# Emission line histograms
let
    SERIES_OPTS = (bins=40.5:0.25:45, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-4.5, :mm), # <-- to reduce wasted space to the left of plot
                   size=(450,550), # <-- controls size of overall plot/image (pixels only)
                   legend = :top)
    accum = Vector{Any}()
    push!(accum, stephist(log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42, label=L"\mathrm{Mg\,II}";  color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Hb  , :Hb_br_norm])        .+ 42, label=L"\mathrm{H\beta}";  color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Ha  , :Ha_br_norm])        .+ 42, label=L"\mathrm{H\alpha}"; color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->"", ylabel=L"\mathrm{Counts}"))
    push!(accum, stephist(log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42, label=L"\mathrm{He\,I}";   color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Pab , :Pab_br_norm])       .+ 42, label=L"\mathrm{Pa}\beta"; color=ith_color(length(accum)+1), SERIES_OPTS..., xlabel=L"\log_{10}(L_{\mathrm{line}}/\mathrm{erg\ s^{-1}})", bottom_margin=(-7., :mm)))
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS...)
    savefig("Figures/Lines_lum.pdf")
end

let
    SERIES_OPTS = (bins=2.5:0.2:5, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-3.5, :mm), # <-- to reduce wasted space to the left of the plot
                   size=(550,650), # <-- controls size of overall plot/image (in pixels only)
                   legend = :topright)
                   
    accum = Vector{Any}()
    push!(accum, stephist(log10.(euclid[qc.MgII, :MgII_2798_br_fwhm]), label=L"\mathrm{Mg\,II}";  color=ith_color(1), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Hb  , :Hb_br_fwhm])       , label=L"\mathrm{H\beta}";  color=ith_color(2), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Ha  , :Ha_br_fwhm])       , label=L"\mathrm{H\alpha}"; color=ith_color(3), SERIES_OPTS..., xformatter=_->"", ylabel=L"\mathrm{Counts}"))
    push!(accum, stephist(log10.(euclid[qc.HeI , :HeI_10832_br_fwhm]), label=L"\mathrm{He\,I}";   color=ith_color(4), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(log10.(euclid[qc.Pab , :Pab_br_fwhm])      , label=L"\mathrm{Pa}\beta"; color=ith_color(5), SERIES_OPTS..., xlabel=L"\log_{10}(\mathrm{FWHM/km\ s^{-1}})", bottom_margin=(-6., :mm)))

    SERIES_OPTS = (SERIES_OPTS..., bins=-500:200.:500)
    push!(accum, stephist(       euclid[qc.MgII, :MgII_2798_br_voff] , label=L"\mathrm{Mg\,II}";  color=ith_color(1), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.Hb  , :Hb_br_voff]        , label=L"\mathrm{H\beta}";  color=ith_color(2), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.Ha  , :Ha_br_voff]        , label=L"\mathrm{H\alpha}"; color=ith_color(3), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.HeI , :HeI_10832_br_voff] , label=L"\mathrm{He\,I}";   color=ith_color(4), SERIES_OPTS..., xformatter=_->""))
    push!(accum, stephist(       euclid[qc.Pab , :Pab_br_voff]       , label=L"\mathrm{Pa}\beta"; color=ith_color(5), SERIES_OPTS..., xlabel=L"\mathrm{V_{off}/km\ s^{-1}}", bottom_margin=(-6., :mm)))

    accum = permutedims(reshape(accum,     div(length(accum), 2), 2)) # reorder subplots so that they appear in the correct order
    plot(reshape(accum, :)..., layout=grid(div(length(accum), 2), 2); GEN_OPTS...)
    savefig("Figures/Lines_FWHM_Voff.pdf")
end


######################################################################
# Bol Lum Mean Histogram
let
    plot(; GEN_OPTS..., legend=:topleft, xlabel=L"\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})", ylabel=L"\mathrm{Counts}")
    vv = filter(!isnan, euclid[qc.good             , :MBH_mean]);  stephist!(vv, label=L"\mathrm{Good\ sample}"; HISTO_OPTS...)
    vv = filter(!isnan, euclid[qc.good .& sub.cosmo, :MBH_mean]);  stephist!(vv, label=L"0.8 < z < 1.9"        ; HISTO_OPTS...)
    savefig("Figures/Mean_BH_Mass.pdf")
end

let
    plot(; GEN_OPTS..., legend=:topleft, xlabel=L"\log_{10}(L_{bol}/\mathrm{erg\ s^{-1}})", ylabel=L"\mathrm{Counts}")
    vv = filter(!isnan, euclid[qc.good             , :Lbol_mean]); stephist!(log10.(vv), label=L"\mathrm{Good\ sample}"; HISTO_OPTS...)
    vv = filter(!isnan, euclid[qc.good .& sub.cosmo, :Lbol_mean]); stephist!(log10.(vv), label=L"0.8 < z < 1.9"        ; HISTO_OPTS...)
    savefig("Figures/Mean_Lbol.pdf")
end


######################################################################
# DESI BH Mass - Euclid Ha BH mass
let
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{BH,\ Euclid,\ H\alpha}} / M_{\mathrm{BH,\ DESI,\ MgII}})", ylabel=L"\mathrm{Counts}")
    ii = qc_match.Ha .& qc_desi.MgII
    vv = filter(!isnan, euclid_match[ii, :MBH_Ha_ShenLiu2012] .- desi[ii, :MBH_MgII_WuShen2022])
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("Figures/BH_mass_cmpDESI_histo.pdf")

    plot(; GEN_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{BH,\ Euclid,\ H\alpha}}/\mathrm{M_{\odot}})", ylabel=L"\log_{10}(M_{\mathrm{BH,\ DESI,\ MgII}}/\mathrm{M_{\odot}})", xlims=(7.5, 10), ylims=(7.5,10))
    ii = qc_match.Ha .& qc_desi.MgII
    scatter!(euclid_match[ii, :MBH_Ha_ShenLiu2012], desi[ii, :MBH_MgII_WuShen2022], label="")
    plot!([xlims()...], [ylims()...], label=L"1:1", linecolor=:black, ls=:dash, lw=2)
    savefig("Figures/BH_mass_cmpDESI_scatter.pdf")
end

let
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(\mathrm{L}_{\mathrm{Euclid,\ H\alpha}}\ /\ \mathrm{L}_{\mathrm{DESI,\ MgII}})", ylabel=L"\mathrm{Counts}")
    ii = qc_match.Ha .& qc_desi.MgII
    vv = log10.(filter(!isnan, euclid_match[ii, :Ha_br_norm] ./ desi[ii, :MgII_2798_br_norm]))
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("Figures/Line_ratio_cmpDESI_Ha.pdf")
end

let
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(\mathrm{L}_{\mathrm{DESI,\ MgII}} / \mathrm{L}_{\mathrm{Euclid,\ H\beta}})", ylabel=L"\mathrm{Counts}")
    ii = qc_match.Hb .& qc_desi.MgII
    vv = log10.(filter(!isnan, desi[ii, :MgII_2798_br_norm] ./ euclid_match[ii, :Hb_br_norm]))
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("Figures/Line_ratio_cmpDESI_Hb.pdf")
end

######################################################################
# Chi2 vs SNR plots
let
    plot(; GEN_OPTS..., xlabel=L"\mathrm{Reduced}\ \chi^2", ylabel=L"\mathrm{Spectrum\ S/N}")
    scatter!(desi[:           , :redchisq], desi[:           , :DER_SNR], label=L"\mathrm{Full\ sample}", ms=2, markerstrokewidth=0)
    scatter!(desi[qc_desi.good, :redchisq], desi[qc_desi.good, :DER_SNR], label=L"\mathrm{Good\ sample}", ms=4, ma=0.6, markershape=:circ)
    plot!(xaxis=:log10, yaxis=:log10)
    savefig("Figures/Chi2vsSNR_cut_DESI.pdf")
end


######################################################################
# SDSS plots
let
    SERIES_OPTS = (ma=0.6,)
    plot(; GEN_OPTS..., xlabel=L"z", ylabel=L"\mathrm{\log_{10}(M_{\mathrm{BH}}/{\mathrm{M_{\odot}}})}", legend=:bottomright)
    scatter!(euclid[qc.Pab , :Redshift], euclid[qc.Pab , :MBH_Pab_Ricci]      , label=L"\mathrm{Pa\beta}", markershape=:dtriangle; SERIES_OPTS...)
    scatter!(euclid[qc.HeI , :Redshift], euclid[qc.HeI , :MBH_HeI_Ricci]      , label=L"\mathrm{He\,I}"  , markershape=:circle   ; SERIES_OPTS...)
    scatter!(euclid[qc.Ha  , :Redshift], euclid[qc.Ha  , :MBH_Ha_ShenLiu2012] , label=L"\mathrm{H\alpha}", markershape=:rect     ; SERIES_OPTS...)
    scatter!(euclid[qc.Hb  , :Redshift], euclid[qc.Hb  , :MBH_Hb_WuShen2022]  , label=L"\mathrm{H\beta}" , markershape=:pentagon ; SERIES_OPTS...)
    scatter!(euclid[qc.MgII, :Redshift], euclid[qc.MgII, :MBH_MgII_WuShen2022], label=L"\mathrm{Mg\,II}" , markershape=:utriangle; SERIES_OPTS...)
    savefig("Figures/MBH_z.png")

    SERIES_OPTS = (ma=0.6,)
    plot(; GEN_OPTS..., legend=:bottomright, xlims=(44, 47.5),
         xlabel=L"\log_{10}(L_{\mathrm{bol}}/\mathrm{erg\ s^{-1}})", ylabel=L"\mathrm{\log_{10}(M_{\mathrm{BH}}/{\mathrm{M_{\odot}}})}")
    scatter!(log10.(euclid[qc.Pab , :Lbol_mean]), euclid[qc.Pab , :MBH_Pab_Ricci]      , label=L"\mathrm{Pa\beta}", markershape=:dtriangle; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.HeI , :Lbol_mean]), euclid[qc.HeI , :MBH_HeI_Ricci]      , label=L"\mathrm{He\,I}"  , markershape=:circle   ; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.Ha  , :Lbol_mean]), euclid[qc.Ha  , :MBH_Ha_ShenLiu2012] , label=L"\mathrm{H\alpha}", markershape=:rect     ; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.Hb  , :Lbol_mean]), euclid[qc.Hb  , :MBH_Hb_WuShen2022]  , label=L"\mathrm{H\beta}" , markershape=:pentagon ; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.MgII, :Lbol_mean]), euclid[qc.MgII, :MBH_MgII_WuShen2022], label=L"\mathrm{Mg\,II}" , markershape=:utriangle; SERIES_OPTS...)
    
    M = 10 .^[ylims()...]
    Ledd = M .* 1.26e38
    for eddratio in [0.1, 1, 10]
        Lbol = eddratio .* Ledd
        plot!(log10.(Lbol), log10.(M), label=L"\eta=%$eddratio", linestyle=:dash)
    end
    plot!()
    savefig("Figures/MBH_Lbol.png")
end


##########################################################################
#Number counts for the paper

@printf("Total number of spectra in Euclid sample: %i \n", nrow(euclid))
@printf("Number of sources in sample with host galaxy contribution: %i \n", count(sub.hostgal))
@printf("Number of sources in sample present in Fu et al.: %i \n", count(sub.FU))
@printf("Number of sources in sample present in DESI: %i \n", count(sub.DESI))
@printf("Number of sources in sample present in QUBRICS: %i \n", count(sub.QUBRICS))
@printf("Good quality spectra in Euclid sample: %i \n", count(qc.good))
@printf("Number of spectra with counterpart in DESI: %i \n", nrow(euclid_match))
@printf("Good quality spectra with counterpart in DESI: %i \n", count(qc_match.good))

println("")

@printf("Mean H magnitude for quality sample: %.1f +/- %.1f \n", mean(euclid[qc.good, :Hmag]), std(euclid[qc.good, :Hmag]))
@printf("Fraction of quality sources with Hmag<21: %.1f%% \n", (length(euclid[qc.good, :Hmag] .< 21)/nrow(euclid))*100)

println("")

iii = findall(sub.FU      .& qc.good)
@printf("Good quality sources present in Fu et al.: %i \n", length(iii))
@printf("	redshift range: %.3f -- %.3f \n", minimum(euclid[iii, :Redshift]), maximum(euclid[iii, :Redshift]))
@printf("	H mag range: %.2f -- %.2f \n", minimum(euclid[iii, :Hmag]),  maximum(euclid[iii, :Hmag]))

println("")

iii = findall(sub.DESI      .& qc.good)
@printf("Good quality sources present in DESI: %i \n", length(iii))
@printf("	redshift range: %.3f -- %.3f \n", minimum(euclid[iii, :Redshift]), maximum(euclid[iii, :Redshift]))
@printf("	H mag range: %.2f -- %.2f \n", minimum(euclid[iii, :Hmag]), maximum(euclid[iii, :Hmag]))

println("")

iii = findall(sub.QUBRICS      .& qc.good)
@printf("Good quality sources present in QUBRICS: %i \n", length(iii))
@printf("	redshift range: %.3f -- %.3f \n", minimum(euclid[iii, :Redshift]), maximum(euclid[iii, :Redshift]))
@printf("	H mag range: %.2f -- %.2f \n", minimum(euclid[iii, :Hmag]), maximum(euclid[iii, :Hmag]))

println("")

@printf("QSO slope for quality sample: %.2f +/- %.2f \n", median(euclid[qc.good, :QSOcont_alpha]), mad(euclid[qc.good, :QSOcont_alpha], normalize=true))
@printf("QSO slope for quality sample at z<0.8: %.2f +/- %.2f \n", median(euclid[qc.good .& sub.lowz, :QSOcont_alpha]), mad(euclid[qc.good .& sub.lowz, :QSOcont_alpha], normalize=true))
@printf("QSO slope for quality sample at 0.8<z<1.9: %.2f +/- %.2f \n", median(euclid[qc.good .& sub.cosmo, :QSOcont_alpha]), mad(euclid[qc.good .& sub.cosmo, :QSOcont_alpha], normalize=true))
@printf("QSO slope for quality sample z>1.9: %.2f +/- %.2f \n", median(euclid[qc.good .& sub.highz, :QSOcont_alpha]), mad(euclid[qc.good .& sub.highz, :QSOcont_alpha], normalize=true))
@printf("QSO slope for DESI quality sample: %.2f +/- %.2f \n", median(desi[qc_desi.good, :QSOcont_alpha]), mad(desi[qc_desi.good, :QSOcont_alpha], normalize=true))

println("")

@printf("Mean luminosity for Ha: %.2f +/- %.2f \n", mean(log10.(euclid[qc.Ha, :Ha_br_norm]) .+ 42), std(log10.(euclid[qc.Ha, :Ha_br_norm]) .+ 42))
@printf("Mean FWHM for Ha: %.2f +/- %.2f \n", mean(log10.(euclid[qc.Ha, :Ha_br_fwhm])), std(log10.(euclid[qc.Ha, :Ha_br_fwhm])))
@printf("Mean v_off for Ha: %.2f +/- %.2f \n", mean(euclid[qc.Ha, :Ha_br_voff]), std(euclid[qc.Ha, :Ha_br_voff]))
@printf("Mean BH mass for Ha: %.2f +/- %.2f \n", mean(euclid[qc.Ha, :MBH_Ha_ShenLiu2012]), std(euclid[qc.Ha, :MBH_Ha_ShenLiu2012]))
@printf("Mean Edd. ratio for Ha: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.Ha, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.Ha, :Edd_ratio]))))

println("")

@printf("Mean luminosity for Hb: %.2f +/- %.2f \n", mean(log10.(euclid[qc.Hb, :Hb_br_norm]) .+ 42), std(log10.(euclid[qc.Hb, :Hb_br_norm]) .+ 42))
@printf("Mean FWHM for Hb: %.2f +/- %.2f \n", mean(log10.(euclid[qc.Hb, :Hb_br_fwhm])), std(log10.(euclid[qc.Hb, :Hb_br_fwhm])))
@printf("Mean v_off for Hb: %.2f +/- %.2f \n", mean(euclid[qc.Hb, :Hb_br_voff]), std(euclid[qc.Hb, :Hb_br_voff]))
@printf("Mean BH mass for Hb: %.2f +/- %.2f \n", mean(filter(!isnan, euclid[qc.Hb, :MBH_Hb_WuShen2022])), std(filter(!isnan, euclid[qc.Hb, :MBH_Hb_WuShen2022])))
@printf("Mean Edd. ratio for Hb: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.Hb, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.Hb, :Edd_ratio]))))

println("")

@printf("Mean luminosity for MgII: %.2f +/- %.2f \n", mean(log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42), std(log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42))
@printf("Mean FWHM for MgII: %.2f +/- %.2f \n", mean(log10.(euclid[qc.MgII, :MgII_2798_br_fwhm])),  std(log10.(euclid[qc.MgII, :MgII_2798_br_fwhm])))
@printf("Mean v_off for MgII: %.2f +/- %.2f \n", mean(euclid[qc.MgII, :MgII_2798_br_voff]), std(euclid[qc.MgII, :MgII_2798_br_voff]))
@printf("Mean BH mass for MgII: %.2f +/- %.2f \n", mean(euclid[qc.MgII, :MBH_MgII_WuShen2022]), std(euclid[qc.MgII, :MBH_MgII_WuShen2022]))
@printf("Mean Edd. ratio for MgII: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.MgII, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.MgII, :Edd_ratio]))))

println("")

@printf("Mean luminosity for HeI: %.2f +/- %.2f \n", mean(log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42), std(log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42))
@printf("Mean FWHM for HeI: %.2f +/- %.2f \n", mean(log10.(euclid[qc.HeI , :HeI_10832_br_fwhm])), std(log10.(euclid[qc.HeI , :HeI_10832_br_fwhm])))
@printf("Mean v_off for HeI: %.2f +/- %.2f \n", mean(euclid[qc.HeI , :HeI_10832_br_voff]), std(euclid[qc.HeI , :HeI_10832_br_voff]))
@printf("Mean BH mass for HeI: %.2f +/- %.2f \n", mean(euclid[qc.HeI, :MBH_HeI_Ricci]), std(euclid[qc.HeI, :MBH_HeI_Ricci]))
@printf("Mean Edd. ratio for HeI: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.HeI, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.HeI, :Edd_ratio]))))

println("")

@printf("Mean luminosity for Pab: %.2f +/- %.2f \n", mean(log10.(euclid[qc.Pab , :Pab_br_norm]) .+ 42), std(log10.(euclid[qc.Pab , :Pab_br_norm]) .+ 42))
@printf("Mean FWHM for Pab: %.2f +/- %.2f \n", mean(log10.(euclid[qc.Pab , :Pab_br_fwhm])), std(log10.(euclid[qc.Pab , :Pab_br_fwhm])))
@printf("Mean v_off for Pab: %.2f +/- %.2f \n", mean(euclid[qc.Pab , :Pab_br_voff]), std(euclid[qc.Pab , :Pab_br_voff]))
@printf("Mean BH mass for Pab: %.2f +/- %.2f \n", mean(euclid[qc.Pab, :MBH_Pab_Ricci]), std(euclid[qc.Pab, :MBH_Pab_Ricci]))
@printf("Mean Edd. ratio for Pab: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.Pab, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.Pab, :Edd_ratio]))))

iii= qc_match.Ha .& qc_desi.MgII
@printf("Mean ratio of Euclid Ha and DESI MgII: %.2f +/- %.2f \n", mean(log10.(euclid_match[iii, :Ha_br_norm] ./ desi[iii, :MgII_2798_br_norm])), std(log10.(euclid_match[iii, :Ha_br_norm] ./ desi[iii, :MgII_2798_br_norm])))

iii = qc_match.Hb .& qc_desi.MgII
@printf("Mean ratio of Euclid Hb and DESI MgII: %.2f +/- %.2f \n", mean(log10.(euclid_match[iii, :Hb_br_norm] ./ desi[iii, :MgII_2798_br_norm])), std(log10.(euclid_match[iii, :Hb_br_norm] ./ desi[iii, :MgII_2798_br_norm])))

println("")

@printf("Mean bolometric luminosity for quality sample: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.good, :Lbol_mean]))), std(filter(!isnan, log10.(euclid[qc.good, :Lbol_mean]))))
@printf("Mean BH mass for quality sample: %.2f +/- %.2f \n", mean(filter(!isnan, euclid[qc.good, :MBH_mean])), std(filter(!isnan, euclid[qc.good, :MBH_mean])))





