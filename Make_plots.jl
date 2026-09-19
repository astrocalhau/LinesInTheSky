using Statistics, StatsBase, Plots, StatsPlots, LaTeXStrings, DataFrames, FITSIO, Printf, QSFit, TypedJSON, SortMerge, GModelFit

######################################################################
# Script to create the various figures for the Euclid Q1 Paper "Lines in the Sky"

mkpath("Figures")
output="Figures"

#f = FITS("results_Euclid_BIC/QSFIT_RESULTS_BIC_w_errors_2.fits")
f = FITS("results_Euclid/QSFIT_RESULTS.fits")
#f = FITS("results_Euclid_earlytype/QSFIT_RESULTS_earlytype.fits")
euclid = DataFrame(f[2])
close(f)

#f = FITS("results_DESI/QSFIT_RESULTS_w_errors.fits")
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
           Ha   = ((df.good .== 1)  .&  (df.Ha_br_reliable        .== 1)), # quality cut for Ha line
           Hb   = ((df.good .== 1)  .&  (df.Hb_br_reliable        .== 1)), # quality cut for Hb line
           MgII = ((df.good .== 1)  .&  (df.MgII_2798_br_reliable .== 1))) # quality cut for MgII line
     
    if "Pab_br_reliable" in names(df)
        out = (out...,
           Pab  = ((df.good .== 1)  .&  (df.Pab_br_reliable       .== 1))) # quality cut for Pab line
    end
    if "HeI_10832_br_reliable" in names(df)
        out = (out...,
           HeI  = ((df.good .== 1)  .&  (df.HeI_10832_br_reliable .== 1))) # quality cut for HeI line
    end
    return out
end
    
qc = quality_cuts(euclid)              # Quality cuts for the whole Euclid catalog
qc_match = quality_cuts(euclid_match)  # Quality cuts for the subset in Euclid catalog having a DESI spectrum
qc_desi = quality_cuts(desi)           # Quality cuts for the DESI catalog

# General plot settings
GEN_OPTS = (fontfamily="Computer Modern", framestyle=:box, grid=false,
            background_color_legend=nothing, foreground_color_legend=nothing, # , legend_column=-1
            #palette = [:red, :green, :blue, :darkorange3, :purple, :gray],
            #palette = [colorant"#648FFF", colorant"#785EF0", colorant"#DC267F", colorant"#FE6100", colorant"#FFB000", :gray],
            palette = [colorant"#a50026", colorant"#fdae61", colorant"#74add1", colorant"#d73027",  colorant"#313695", :gray],
            thickness_scaling=1.4)
# , titlefontsize=8, guidefontsize=8, tickfontsize=8, legendfontsize=7

HISTO_OPTS = (alpha=0, fillalpha=0.8, fill=true, linewidth=false)
 
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
    SERIES_OPTS = (bins=0:0.25:maximum(euclid.Redshift), fillalpha=0.7)
    stephist!(euclid[:          , :Redshift], label=L"\mathrm{Full\ sample}"            , fill=false, color=:black; SERIES_OPTS...)
    stephist!(euclid[sub.FU     , :Redshift], label=L"\mathrm{Fu\ et\ al.}\ (2026)", fill= true, linewidth=false, alpha=0, color=ith_color(1); SERIES_OPTS...)
    stephist!(euclid[sub.DESI   , :Redshift], label=L"\mathrm{DESI}"                    , fill=true , linewidth=false, alpha=0, color=ith_color(2);  SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS, :Redshift], label=L"\mathrm{QUBRICS}"                 , fill=true , linewidth=false, alpha=0, color=ith_color(3), bins=0:0.25:maximum(euclid.Redshift), fillalpha=1)
    savefig("Figures/Redshift_Hist_full_sample.pdf")

    plot(; GEN_OPTS..., title="", xlabel=L"z", ylabel=L"\mathrm{Counts}")
    stephist!(euclid[                 qc.good, :Redshift], label=L"\mathrm{Good\ sample}"             , fill=false, color=:black; SERIES_OPTS...)
    stephist!(euclid[sub.FU       .&  qc.good, :Redshift], label=L"\mathrm{Fu\ et\ al.}\ (2026)" , fill= true, linewidth=false, alpha=0, color=ith_color(1); SERIES_OPTS...)
    stephist!(euclid[sub.DESI     .&  qc.good, :Redshift], label=L"\mathrm{DESI}"                     , fill= true, linewidth=false, alpha=0, color=ith_color(2); SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS  .&  qc.good, :Redshift], label=L"\mathrm{QUBRICS}"                  , fill= true, linewidth=false, alpha=0, color=ith_color(3); SERIES_OPTS...)
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
    vv = filter(!isnan, euclid[qc.MgII, :MBH_MgII_WuShen2022]); push!(accum, stephist(vv, label=""   , color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Mg\,II}", 9, :dark, rotation=0))
    vv = filter(!isnan, euclid[qc.Hb  , :MBH_Hb_WuShen2022])  ; push!(accum, stephist(vv, label="" , color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\beta}", 9, :dark, rotation=0))
    vv = filter(!isnan, euclid[qc.Ha  , :MBH_Ha_ShenLiu2012]) ; push!(accum, stephist(vv, label="", color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"", ylabel=L"\mathrm{Counts}")); add_μσ!(vv); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\alpha}", 9, :dark, rotation=0))
    vv = filter(!isnan, euclid[qc.HeI , :MBH_HeI_Ricci])      ; push!(accum, stephist(vv, label=""  , color=ith_color(length(accum)+1); SERIES_OPTS..., xformatter=_->"")); add_μσ!(vv); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{He\,I}", 9, :dark, rotation=0))
    vv = filter(!isnan, euclid[qc.Pab , :MBH_Pab_Ricci])      ; push!(accum, stephist(vv, label="", color=ith_color(length(accum)+1); SERIES_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{BH}}/{M_{\odot}})", bottom_margin=(-1., :mm))); add_μσ!(vv); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Pa}\beta", 9, :dark, rotation=0))
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS..., legend=:topright)
    savefig("$(output)/BH_mass.pdf")

    # Comparison between Ha and Hb
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{H}\alpha}/M_{\mathrm{H}\beta})", ylabel=L"\mathrm{Counts}")
    ii = findall(qc.Ha  .&  qc.Hb)
    vv = filter(!isnan, euclid[ii, :MBH_Ha_ShenLiu2012] .- euclid[ii, :MBH_Hb_WuShen2022])
    stephist!(vv, label="", bins=minimum(vv):0.25:maximum(vv); HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("$(output)/BH_mass_Ha_vs_Hb.pdf")
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
            out = stephist(vv, label="", color=color; SERIES_OPTS..., kws...)
            vline!([median(vv)], label="", color=:black, linewidth=2, linestyle=:dash)
            annotate!((0.7, 0.8), fontfamily="Computer Modern", Plots.text(label *" "* L"\ (\tilde{\mu}=%$sμ)", 9, :dark, rotation=0))
        else
            out = stephist(vv, label=label, color=color; SERIES_OPTS..., kws...)
        end
    end
    
    # plot(; GEN_OPTS..., legend=:topright, xlabel=L"\alpha_{\lambda}", ylabel=L"\mathrm{Counts}")
    SERIES_OPTS = (bins=minimum(euclid.QSOcont_alpha):0.25:maximum(euclid.QSOcont_alpha), fillalpha=0.8, fill=true, linewidth=false,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-2.5, :mm), # <-- to reduce wasted space to the left of the plot
                   size=(500, 850)) # <-- Controls size of overall plot
    accum = Vector{Any}()
    push!(accum, add_series!(euclid[qc.good              , :QSOcont_alpha], L"\mathrm{Good\ sample}"  , :black      , SERIES_OPTS, xformatter=_->"", fill=false, linewidth=4))
    push!(accum, add_series!(euclid[qc.good .&  sub.lowz , :QSOcont_alpha], L"z < 0.8"                , ith_color(1), SERIES_OPTS, alpha =0, xformatter=_->""))
    push!(accum, add_series!(euclid[qc.good .&  sub.cosmo, :QSOcont_alpha], L"0.8 < z < 1.9"          , ith_color(2), SERIES_OPTS, alpha =0, xformatter=_->"", ylabel=L"\mathrm{Counts}"))
    push!(accum, add_series!(euclid[qc.good .&  sub.highz, :QSOcont_alpha], L"z > 1.9"                , ith_color(3), SERIES_OPTS, alpha =0, xformatter=_->""))
    push!(accum, add_series!(desi[qc_desi.good           , :QSOcont_alpha], L"\mathrm{DESI\,\,spectra}"                   , ith_color(4), SERIES_OPTS,alpha =0, xlabel=L"\alpha_{\lambda}", bottom_margin=(-6., :mm)))
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS..., legend=:topright)
    savefig("$(output)/QSOcont_alpha_Hist_Redshift.pdf")
end


######################################################################
# Chi2 vs SNR plots
let
    plot(; GEN_OPTS..., xlabel=L"\mathrm{Reduced}\ \chi^2", ylabel=L"\mathrm{Spectrum\ S/N}")
    scatter!(euclid[:      , :redchisq], euclid[:      , :DER_SNR], label=L"\mathrm{Full\ sample}", ms=2, markerstrokewidth=0)
    scatter!(euclid[qc.good, :redchisq], euclid[qc.good, :DER_SNR], label=L"\mathrm{Good\ sample}", ms=4, ma=0.6, markershape=:circ)
    plot!(xaxis=:log10, yaxis=:log10)
    savefig("$(output)/Chi2vsSNR_cut.pdf")
end


######################################################################
# H magnitude histogram
let
    SERIES_OPTS = (bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)),  alpha=1, fillalpha=0.7)
    plot(; GEN_OPTS..., xlabel=L"H_E\ \mathrm{magnitude}", ylabel=L"\mathrm{Counts}", title=L"\mathrm{\textbf{Good\ sample}}", legend=:topleft)
    stephist!(euclid[               qc.good, :Hmag], label=L"\mathrm{Merged\ sample}"           , fill=false, color=:black; SERIES_OPTS...)
    stephist!(euclid[sub.FU      .& qc.good, :Hmag], label=L"\mathrm{Fu\ et\ al.\ (in\ prep.)}" , fill= true, color=ith_color(1); SERIES_OPTS...)
    stephist!(euclid[sub.DESI    .& qc.good, :Hmag], label=L"\mathrm{DESI}"                     , fill= true, color=ith_color(2); SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS .& qc.good, :Hmag], label=L"\mathrm{QUBRICS}"                  , fill= true, color=ith_color(3); SERIES_OPTS...)
    stephist!(euclid[               qc.good, :Hmag], label=""                                   , fill=false, color=:black; SERIES_OPTS...)
    savefig("Figures/Hmag_Hist.pdf")
end

# Hmag hist vs scatter vs z hist
let
    SCATTER_OPTS = (ms=4, ma=1.0)
    HIST_OPTS = (bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)), alpha=1, fillalpha=0.7, linewidth=false)
    SERIES_OPTS = (alpha=1, fillalpha=0., linewidth=2)
   
            
    histo = stephist(euclid[                 qc.good, :Redshift], label=L"\mathrm{Good\ sample}"             , fill=false, color=:black,  yticks=([0, 100, 200, 300],[L"0", L"100", L"200", L"300"]), xformatter=Returns(""), grid=false, framestyle=:box, legendfontsize=8, foreground_color_legend = nothing, background_color_legend = nothing; SERIES_OPTS...)
    stephist!(euclid[sub.FU       .&  qc.good, :Redshift],        label=""             , fill= true, color=ith_color(1); SERIES_OPTS...)
    stephist!(euclid[sub.DESI     .&  qc.good, :Redshift],        label=""             , fill= true, color=ith_color(2); SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS  .&  qc.good, :Redshift],        label=""             , fill= true, color=ith_color(3); SERIES_OPTS...)
    
    dots = (i = findall(sub.FU      .& qc.good); scatter(euclid[sub.FU      .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{Fu\ et\ al.}\ (2026)", markershape=:rect, xlabel=L"z", ylabel=L"H_\mathrm{E}", ylims=(14, 23), formatter=:latex, legend=:bottomright; SCATTER_OPTS..., GEN_OPTS...))
    i = findall(sub.DESI    .& qc.good); scatter!(euclid[sub.DESI    .& qc.good, :Redshift], euclid[i, :Hmag],  label=L"\mathrm{DESI}"       , markershape=:pentagon ; SCATTER_OPTS...)
    i = findall(sub.QUBRICS .& qc.good); scatter!(euclid[sub.QUBRICS .& qc.good, :Redshift], euclid[i, :Hmag],  label=L"\mathrm{QUBRICS}"    , markershape=:utriangle; SCATTER_OPTS...)
    
    histo2 = stephist(euclid[               qc.good, :Hmag], label=""           , fill=false, color=:black, permute = (:x, :y), xlims=(14, 23), yticks=([0, 125, 250],[L"0", L"125", L"250"]), bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)), xformatter=Returns(""), grid=false, framestyle=:box, margin = -10Plots.px; SERIES_OPTS...)
    stephist!(euclid[sub.FU      .& qc.good, :Hmag], label="" , fill= true, color=ith_color(1), permute = (:x, :y), xlims=(14, 23), bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)); SERIES_OPTS...)
    stephist!(euclid[sub.DESI    .& qc.good, :Hmag], label=""                     , fill= true, color=ith_color(2), permute = (:x, :y), xlims=(14, 23), bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)); SERIES_OPTS...)
    stephist!(euclid[sub.QUBRICS .& qc.good, :Hmag], label=""                  , fill= true, color=ith_color(3), permute = (:x, :y), xlims=(14, 23), bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)); SERIES_OPTS...)
    stephist!(euclid[               qc.good, :Hmag], label=""                                   , fill=false, color=:black, permute = (:x, :y), xlims=(14, 23), bins=minimum(filter(!isnan, euclid.Hmag)):0.25:maximum(filter(!isnan, euclid.Hmag)); SERIES_OPTS...)
    
    plot(histo, plot(framestyle=:none), dots, histo2, link=:both, layout=grid(2,2, heights=(1/4,3/4), widths=(5/6, 1/6)), thickness_scaling=1.4, top_margin=[-5Plots.mm 0Plots.mm], bottom_margin=[0Plots.mm -5Plots.mm], left_margin=[-2Plots.mm -10Plots.mm])
    
    savefig("Figures/Hmag_redshift_scatter_w_hitst.pdf")
end


######################################################################
# Hmag vs redshift
let
    SERIES_OPTS = (ms=6, ma=0.6, legend=:bottomright)
    plot(; GEN_OPTS..., title=L"\mathrm{\textbf{Good\ sample}}", xlabel=L"z", ylabel=L"H_\mathrm{E}", ylims=(14, 23))
    i = findall(sub.FU      .& qc.good); scatter!(euclid[sub.FU      .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{Fu\ et\ al.\ (2026)}", markershape=:rect     ; SERIES_OPTS...)
    i = findall(sub.DESI    .& qc.good); scatter!(euclid[sub.DESI    .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{DESI}"       , markershape=:pentagon ; SERIES_OPTS...)
    i = findall(sub.QUBRICS .& qc.good); scatter!(euclid[sub.QUBRICS .& qc.good, :Redshift], euclid[i, :Hmag], label=L"\mathrm{QUBRICS}"    , markershape=:utriangle; SERIES_OPTS...)
    savefig("$(output)/Hmag_vs_z_cut.pdf")
end


######################################################################
# Emission line histograms
let
    
    function add_μσ_lum!(vv, color; x=0.03, y1=0.8, y2=0.6)
        μ = median(vv)
        σ = mad(vv, normalize=true)
        sμ = @sprintf("%.2f", μ)
        sσ = @sprintf("%.2f", σ)
        vline!([median(vv)]              , label="", color=:black, linestyle=:dash, linewidth=3)
    end
    
    SERIES_OPTS = (bins=40.5:0.25:45, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-4.5, :mm), # <-- to reduce wasted space to the left of plot
                   size=(450,550), # <-- controls size of overall plot/image (pixels only)
                   legend = :top)
    accum = Vector{Any}()
    vv = log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42; push!(accum, stephist(vv, label="";  color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(length(accum))); annotate!((0.5, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Mg\,II}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.Hb  , :Hb_br_norm])        .+ 42; push!(accum, stephist(vv, label="";  color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(length(accum))); annotate!((0.5, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\beta}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.Ha  , :Ha_br_norm])        .+ 42; push!(accum, stephist(vv, label=""; color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->"", ylabel=L"\mathrm{Counts}"));add_μσ_lum!(vv, ith_color(length(accum))); annotate!((0.5, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\alpha}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42; push!(accum, stephist(vv, label="";   color=ith_color(length(accum)+1), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(length(accum)), x=0.7); annotate!((0.5, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{He\,I}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.Pab , :Pab_br_norm])       .+ 42; push!(accum, stephist(vv, label=""; color=ith_color(length(accum)+1), SERIES_OPTS..., xlabel=L"\log_{10}(L_{\mathrm{line}}/\mathrm{erg\ s^{-1}})", bottom_margin=(-7., :mm)));add_μσ_lum!(vv, ith_color(length(accum)), x=0.7); annotate!((0.5, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Pa}\beta", 9, :dark, rotation=0))
    plot(accum..., layout=grid(length(accum), 1); GEN_OPTS...)
    savefig("$(output)/Lines_lum.pdf")
end

let

    function add_μσ_lum!(vv, color; x=0.03, y1=0.8, y2=0.6)
        μ = median(vv)
        σ = mad(vv, normalize=true)
        sμ = @sprintf("%.2f", μ)
        sσ = @sprintf("%.2f", σ)
        vline!([median(vv)]              , label="", color=:black, linestyle=:dash, linewidth=3)
    end
    
    SERIES_OPTS = (bins=2.5:0.2:5, HISTO_OPTS...,
                   bottom_margin=(-3.5, :mm), top_margin=(-1.5, :mm), # <-- these are necessary to reduce space between the subplots
                   left_margin=(-3.5, :mm), # <-- to reduce wasted space to the left of the plot
                   size=(550,650)) # <-- controls size of overall plot/image (in pixels only)
                   
                   
    accum = Vector{Any}()
    vv = log10.(euclid[qc.MgII, :MgII_2798_br_fwhm]); push!(accum, stephist(vv, label="";  color=ith_color(1), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(1)); annotate!((0.85, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Mg\,II}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.Hb  , :Hb_br_fwhm]) ; push!(accum, stephist(vv       , label="";  color=ith_color(2), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(2)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\beta}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.Ha  , :Ha_br_fwhm]) ; push!(accum, stephist(vv       , label=""; color=ith_color(3), SERIES_OPTS..., xformatter=_->"", ylabel=L"\mathrm{Counts}"));add_μσ_lum!(vv, ith_color(3)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\alpha}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.HeI , :HeI_10832_br_fwhm]); push!(accum, stephist(vv, label="";   color=ith_color(4), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(4)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{He\,I}", 9, :dark, rotation=0))
    vv = log10.(euclid[qc.Pab , :Pab_br_fwhm]); push!(accum, stephist(vv      , label=""; color=ith_color(5), SERIES_OPTS..., xlabel=L"\log_{10}(\mathrm{FWHM/km\ s^{-1}})", bottom_margin=(-6., :mm)));add_μσ_lum!(vv, ith_color(5)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Pa}\beta", 9, :dark, rotation=0))

    SERIES_OPTS = (SERIES_OPTS..., bins=-500:200.:500)
    vv = euclid[qc.MgII, :MgII_2798_br_voff]; push!(accum, stephist(       vv , label="";  color=ith_color(1), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(1)); annotate!((0.85, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Mg\,II}", 9, :dark, rotation=0))
    vv = euclid[qc.Hb  , :Hb_br_voff]; push!(accum, stephist(       vv        , label="";  color=ith_color(2), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(2)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\beta}", 9, :dark, rotation=0))
    vv = euclid[qc.Ha  , :Ha_br_voff]; push!(accum, stephist(       vv        , label=""; color=ith_color(3), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(3)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{H\alpha}", 9, :dark, rotation=0))
    vv = euclid[qc.HeI , :HeI_10832_br_voff]; push!(accum, stephist(       vv , label="";   color=ith_color(4), SERIES_OPTS..., xformatter=_->""));add_μσ_lum!(vv, ith_color(4)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{He\,I}", 9, :dark, rotation=0))
    vv = euclid[qc.Pab , :Pab_br_voff]; push!(accum, stephist(       vv       , label=""; color=ith_color(5), SERIES_OPTS..., xlabel=L"v_\mathrm{off} \mathrm{/km\ s^{-1}}", bottom_margin=(-6., :mm)));add_μσ_lum!(vv, ith_color(5)); annotate!((0.9, 0.8), fontfamily="Computer Modern", Plots.text(L"\mathrm{Pa}\beta", 9, :dark, rotation=0))

    accum = permutedims(reshape(accum,     div(length(accum), 2), 2)) # reorder subplots so that they appear in the correct order
    plot(reshape(accum, :)..., layout=grid(div(length(accum), 2), 2); GEN_OPTS...)
    savefig("$(output)/Lines_FWHM_Voff.pdf")
end


######################################################################
# Bol Lum Mean Histogram
let
    plot(; GEN_OPTS..., legend=:topleft, xlabel=L"\log_{10}(M_{\mathrm{BH}}/{M_{\odot}})", ylabel=L"\mathrm{Counts}")
    vv = filter(!isnan, euclid[qc.good             , :MBH_mean]);  stephist!(vv, label=L"\mathrm{Good\ sample}"; HISTO_OPTS...)
    vv = filter(!isnan, euclid[qc.good .& sub.cosmo, :MBH_mean]);  stephist!(vv, label=L"0.8 < z < 1.9"        ; HISTO_OPTS...)
    savefig("$(output)/Mean_BH_Mass.pdf")
end

let
    plot(; GEN_OPTS..., legend=:topleft, xlabel=L"\log_{10}(L_{\mathrm{bol}}/\mathrm{erg\ s^{-1}})", ylabel=L"\mathrm{Counts}")
    vv = filter(!isnan, euclid[qc.good             , :Lbol_mean]); stephist!(log10.(vv), label=L"\mathrm{Good\ sample}"; HISTO_OPTS...)
    vv = filter(!isnan, euclid[qc.good .& sub.cosmo, :Lbol_mean]); stephist!(log10.(vv), label=L"0.8 < z < 1.9"        ; HISTO_OPTS...)
    savefig("$(output)/Mean_Lbol.pdf")
end


######################################################################
# DESI BH Mass - Euclid Ha BH mass
let
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{BH},\ Euclid,\ \mathrm{H\alpha}}/M_{\mathrm{BH,\ DESI,\ MgII}})", ylabel=L"\mathrm{Counts}")
    ii = qc_match.Ha .& qc_desi.MgII
    vv = filter(!isnan, euclid_match[ii, :MBH_Ha_ShenLiu2012] .- desi[ii, :MBH_MgII_WuShen2022])
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("$(output)/BH_mass_cmpDESI_histo.pdf")

    plot(; GEN_OPTS..., xlabel=L"\log_{10}(M_{\mathrm{BH,}\ Euclid,\ \mathrm{H\alpha}}/{M_{\odot}})", ylabel=L"\log_{10}(M_{\mathrm{BH,\ DESI,\ MgII}}/{M_{\odot}})",xlims=(6.8, 10.1), ylims=(7,10.1), legend=:topleft, legendfontsize=8)
    ii = qc_match.Ha .& qc_desi.MgII
    iii = (abs.(euclid_match.MBH_Ha_ShenLiu2012 .- desi.MBH_MgII_WuShen2022) .> 0.5) .& qc_match.Ha .& qc_desi.MgII
    iv = ((abs.(euclid_match.Ha_br_center .- (12500 ./(1 .+ euclid_match.Redshift))) .< 150) .| (abs.(euclid_match.Ha_br_center .- (18500 ./(1 .+ euclid_match.Redshift))) .< 150)) .& (abs.(euclid_match.MBH_Ha_ShenLiu2012 .- desi.MBH_MgII_WuShen2022) .> 0.5) .& qc_match.Ha .& qc_desi.MgII
    
    
    #scatter!(euclid_match[ii, :MBH_Ha_ShenLiu2012], desi[ii, :MBH_MgII_WuShen2022], yerr=desi[ii, :MBH_MgII_WuShen2022_error], xerr=euclid_match[ii, :MBH_Ha_ShenLiu2012_error], label="")
    scatter!(euclid_match[ii, :MBH_Ha_ShenLiu2012], desi[ii, :MBH_MgII_WuShen2022], label="")
    scatter!(euclid_match[iii, :MBH_Ha_ShenLiu2012], desi[iii, :MBH_MgII_WuShen2022], label=L"|\mathrm{log_{10}}(M_{\mathrm{BH,}\,Euclid})-\mathrm{log_{10}}(M_{\mathrm{BH, \, DESI}})|>0.5")
    scatter!(euclid_match[iv, :MBH_Ha_ShenLiu2012], desi[iv, :MBH_MgII_WuShen2022], label=L"\mathrm{Line \,\, centre} - \mathrm{spectrum \,\, edge}<150 \AA", markershape=:dtriangle)
    #plot!([9.5], [7.5], yerr=[0.05], xerr=[0.24], label="")
    plot!([xlims()...], [ylims()...], label=L"1:1", linecolor=:black, ls=:dash, lw=2)
    savefig("$(output)/BH_mass_cmpDESI_scatter.pdf")
end

#######################################################################
# Line ratios
# Within Euclid
let
    plot(; GEN_OPTS..., xlabel=L"L_{\mathrm{H\alpha}}\ /\ L_{\mathrm{H\beta}}", ylabel=L"\mathrm{Counts}")
    ii = qc.Ha .& qc.Hb
    vv = (filter(!isnan, euclid[ii, :Ha_br_norm] ./ euclid[ii, :Hb_br_norm]))
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("$(output)/Line_ratio_Ha_Hb.pdf")
end

# Comparison with DESI
let
    plot(; GEN_OPTS..., xlabel=L"\log_{10}(L_{\mathrm{Euclid,\ H\alpha}}\ /\ L_{\mathrm{DESI,\ MgII}})", ylabel=L"\mathrm{Counts}")
    ii = qc_match.Ha .& qc_desi.MgII
    vv = log10.(filter(!isnan, euclid_match[ii, :Ha_br_norm] ./ desi[ii, :MgII_2798_br_norm]))
    stephist!(vv, label=""; HISTO_OPTS...)
    add_μσ!(vv, y1=0.9, y2=0.8)
    savefig("$(output)/Line_ratio_cmpDESI_Ha.pdf")
end

let
    plot(; GEN_OPTS..., xlabel=L"L_{\mathrm{DESI,\ MgII}} / L_{\mathrm{Euclid,\ H\beta}}", ylabel=L"\mathrm{Counts}")
    ii = qc_match.Hb .& qc_desi.MgII
    vv = (filter(!isnan, desi[ii, :MgII_2798_br_norm] ./ euclid_match[ii, :Hb_br_norm]))
    stephist!(vv, bins=0:0.55:3.5, label=""; HISTO_OPTS...)
    add_μσ!(vv, x=0.8, y1=0.9, y2=0.8)
    savefig("$(output)/Line_ratio_cmpDESI_Hb.pdf")
end

######################################################################
# Chi2 vs SNR plots
let
    plot(; GEN_OPTS..., xlabel=L"\mathrm{Reduced}\ \chi^2", ylabel=L"\mathrm{Spectrum\ S/N}")
    scatter!(desi[:           , :redchisq], desi[:           , :DER_SNR], label=L"\mathrm{Full\ sample}", ms=2, markerstrokewidth=0)
    scatter!(desi[qc_desi.good, :redchisq], desi[qc_desi.good, :DER_SNR], label=L"\mathrm{Good\ sample}", ms=4, ma=0.6, markershape=:circ)
    plot!(xaxis=:log10, yaxis=:log10)
    savefig("$(output)/Chi2vsSNR_cut_DESI.pdf")
end


######################################################################
# SDSS plots
let
    SERIES_OPTS = (ma=0.6,)
    plot(; GEN_OPTS..., xlabel=L"z", ylabel=L"\mathrm{\log_{10}(M_{\mathrm{BH}}/{M_{\odot}})}", legend=:bottomright)
    scatter!(euclid[qc.Pab , :Redshift], euclid[qc.Pab , :MBH_Pab_Ricci]      , label=L"\mathrm{Pa\beta}", markershape=:dtriangle; SERIES_OPTS...)
    scatter!(euclid[qc.HeI , :Redshift], euclid[qc.HeI , :MBH_HeI_Ricci]      , label=L"\mathrm{He\,I}"  , markershape=:circle   ; SERIES_OPTS...)
    scatter!(euclid[qc.Ha  , :Redshift], euclid[qc.Ha  , :MBH_Ha_ShenLiu2012] , label=L"\mathrm{H\alpha}", markershape=:rect     ; SERIES_OPTS...)
    scatter!(euclid[qc.Hb  , :Redshift], euclid[qc.Hb  , :MBH_Hb_WuShen2022]  , label=L"\mathrm{H\beta}" , markershape=:pentagon ; SERIES_OPTS...)
    scatter!(euclid[qc.MgII, :Redshift], euclid[qc.MgII, :MBH_MgII_WuShen2022], label=L"\mathrm{Mg\,II}" , markershape=:utriangle; SERIES_OPTS...)
    savefig("$(output)/MBH_z.png")

    SERIES_OPTS = (ma=0.6,)
    plot(; GEN_OPTS..., legend=:bottomright, xlims=(44, 47.5),
         xlabel=L"\log_{10}(L_{\mathrm{bol}}/\mathrm{erg\ s^{-1}})", ylabel=L"\mathrm{\log_{10}(M_{\mathrm{BH}}/{M_{\odot}})}")
    scatter!(log10.(euclid[qc.Pab , :Lbol_mean]), euclid[qc.Pab , :MBH_Pab_Ricci]      , label=L"\mathrm{Pa\beta}", markershape=:dtriangle; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.HeI , :Lbol_mean]), euclid[qc.HeI , :MBH_HeI_Ricci]      , label=L"\mathrm{He\,I}"  , markershape=:circle   ; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.Ha  , :Lbol_mean]), euclid[qc.Ha  , :MBH_Ha_ShenLiu2012] , label=L"\mathrm{H\alpha}", markershape=:rect     ; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.Hb  , :Lbol_mean]), euclid[qc.Hb  , :MBH_Hb_WuShen2022]  , label=L"\mathrm{H\beta}" , markershape=:pentagon ; SERIES_OPTS...)
    scatter!(log10.(euclid[qc.MgII, :Lbol_mean]), euclid[qc.MgII, :MBH_MgII_WuShen2022], label=L"\mathrm{Mg\,II}" , markershape=:utriangle; SERIES_OPTS...)
    
    M = 10 .^[ylims()...]
    Ledd = M .* 1.26e38
    for eddratio in [0.1, 1, 10]
        Lbol = eddratio .* Ledd
        plot!(log10.(Lbol), log10.(M), label=L"\lambda=%$eddratio", linestyle=:dash)
    end
    plot!()
    savefig("$(output)/MBH_Lbol.png")


    SERIES_OPTS = (ma=0.6,)
    plot(; GEN_OPTS..., legend=:bottomright, ylims=(44, 47.5), xlims=(6,10),
         ylabel=L"\log_{10}(L_{\mathrm{bol}}/\mathrm{erg\ s^{-1}})", xlabel=L"\mathrm{\log_{10}(M_{\mathrm{BH}}/{M_{\odot}})}")
    scatter!(euclid[qc.Pab , :MBH_Pab_Ricci], log10.(euclid[qc.Pab , :Lbol_mean])      , label=L"\mathrm{Pa\beta}", markershape=:dtriangle; SERIES_OPTS...)
    scatter!(euclid[qc.HeI , :MBH_HeI_Ricci], log10.(euclid[qc.HeI , :Lbol_mean])      , label=L"\mathrm{He\,I}"  , markershape=:circle   ; SERIES_OPTS...)
    scatter!(euclid[qc.Ha  , :MBH_Ha_ShenLiu2012], log10.(euclid[qc.Ha  , :Lbol_mean]) , label=L"\mathrm{H\alpha}", markershape=:rect     ; SERIES_OPTS...)
    scatter!(euclid[qc.Hb  , :MBH_Hb_WuShen2022], log10.(euclid[qc.Hb  , :Lbol_mean])  , label=L"\mathrm{H\beta}" , markershape=:pentagon ; SERIES_OPTS...)
    scatter!(euclid[qc.MgII, :MBH_MgII_WuShen2022], log10.(euclid[qc.MgII, :Lbol_mean]), label=L"\mathrm{Mg\,II}" , markershape=:utriangle; SERIES_OPTS...)
    
    M = 10 .^[xlims()...]
    Ledd = M .* 1.26e38
    for eddratio in [0.1, 1, 10]
        Lbol = eddratio .* Ledd
        plot!(log10.(M), log10.(Lbol), label=L"\lambda=%$eddratio", linestyle=:dash)
    end
    plot!(legend=:topleft)
    savefig("$(output)/Lbol_MBH.png")


end

#########################################################################
#########################################################################
#Ploting spectra 
#Fig. A.1

let
    SERIES_OPTS = (size=(1150,350), GEN_OPTS...,
    		bottom_margin=(3, :mm),
    		xlims=(1200.0,1850.0))
    
    #First panel    
    spectra = FITS("input_Euclid/fits/-602351047484917583.fits")
    spec = DataFrame(spectra[2])
    close(spectra)
    
    mask = (spec.WAVELENGTH .> 12000) .& (spec.WAVELENGTH .< 18500) #to avoid edges of spectra, which are always problematic
    
    SPE = plot(spec.WAVELENGTH[mask] ./10               , spec.SIGNAL[mask]             , title=L"-602351047484917583"               , label=""           , linecolor =:darkblue,     linestyle=:solid    , ylims=(0, 1.5))
    
    #lines
    Ha =     vline!([1492.0]                    , label=""                , color="black"         , linewidth=1        , linestyle=:dash)
    OIII1 =  vline!([1525.0]                    , label=""                , color="black"         , linewidth=1        , linestyle=:dash)
    OIII2 =  vline!([1538.0]                    , label=""                , color="black"         , linewidth=1        , linestyle=:dash)
    Hgamma = vline!([1332.0]                    , label=""                , color="black"         , linewidth=1        , linestyle=:dash)
    
    #Annotations for the lines
    Ha_label = annotate!(1480.0                 , 1.2                     , Plots.text(L"\rm H\beta"                       , 8                         , :dark            , rotation=90))
    OIII1_label = annotate!(1509.0              , 1.23                    , Plots.text(L"${\rm O\textsc{iii}\lambda 4959}$", 8                         , :dark            , rotation=90))
    OIII2_label = annotate!(1553.0              , 1.23                    , Plots.text(L"${\rm O\textsc{iii}\lambda 5007}$", 8                         , :dark            , rotation=90))
    Hgamma_label = annotate!(1310.0             , 1.2                     , Plots.text(L"\rm H\gamma"                      , 8                         , :dark            , rotation=90))
    
    

    #Second panel
    spectra = FITS("input_Euclid/fits/-623046191470992718.fits")
    spec = DataFrame(spectra[2])
    close(spectra)
    
    SPE2 = plot(spec.WAVELENGTH ./10               , spec.SIGNAL             , title=L"-623046191470992718"               , label=""                     , linecolor = :darkblue,     linestyle=:solid    , ylims=(0, 0.2))
    
    #Assmbling the panels together
    plot(SPE, SPE2, layout=grid(1, 2); SERIES_OPTS..., xlabel=L"\mathrm{Wavelength\,\, [nm]}", ylabel=L"\mathrm{Flux \, \, [10^{-16} \,\, erg \,\, cm^{-2} \,\, s^{-1}]}")
    savefig("$(output)/Example_Spectrum.pdf")
end 

#########################################################################
#Fig. A.2

let
    SERIES_OPTS = (size=(1150,350), GEN_OPTS...,
    		bottom_margin=(3, :mm),
    		xlims=(1200,1850))
    
    accum = Vector{Any}()
    filename=["2689906553676437265", "2721144714662952585"]
    
    for file in filename
        spectra = FITS("input_Euclid/fits/$(file).fits")
        spec = DataFrame(spectra[2])
        close(spectra)
        
        mask = (spec.WAVELENGTH .> 12000) .& (spec.WAVELENGTH .< 18500) #to avoid edges of spectra, which are always problematic
        
        push!(accum, plot(spec.WAVELENGTH[mask] ./10             , spec.SIGNAL[mask]             , title=latexstring("$(file)")       , label=""        , linecolor = :darkblue, linestyle=:solid))
    end
    
    plot(accum..., layout=grid(1, length(accum)); SERIES_OPTS..., xlabel=L"\mathrm{Wavelength\,\, [nm]}", ylabel=L"\mathrm{Flux \, \, [10^{-16} \,\, erg \,\, cm^{-2} \,\, s^{-1}]}")
    savefig("$(output)/Bad_spectra_example.pdf")
end 

##########################################################################
# Fig. A.1
let
    SERIES_OPTS = (size=(1200,1600), GEN_OPTS...,
    		bottom_margin=(1, :mm),
    		top_margin=(1,:mm))
    
    accum = Vector{Any}()
    filename = [-663458323479296921, -640749119461755401, -646101801481095171, -603446588505854734]
    
    for file in filename
        ######
        # Read the JSON files
        res = TypedJSON.deserialize("results_Euclid/JSON/$(file).json.gz")
               
        waveor = GModelFit.coords(domain(res.data)) ./10 #original wavelength input (restframed)
        wavecomp = GModelFit.coords(domain(res.bestfit)) ./10 #Domain for the components
        fluxor = values(res.data) #original flux input (restframed)
        flux = GModelFit.folded(res.bestfit)
        broad = res.bestfit(:BroadLines) #Broad lines component
        iron = res.bestfit(:Iron) #Iron lines component
        balmer = res.bestfit(:Balmer)
        #galax = res.bestfit(:Galaxy)
        QSOcont = res.bestfit(:QSOcont) #QSO continuum component
        Nuisance = res.bestfit(:NuisanceLines) #Nuisance lines component

        esp = plot(waveor          , fluxor            , label=L"{\mathrm{Spectrum}}"                        , linestyle=:solid             , linecolor=:black                          )
        plot!(waveor               , flux              , label=L"\mathrm{Model}"                             , linestyle=:solid             , linecolor=ith_color(1)             , linewidth=2)
        plot!(wavecomp             , broad             , label=L"\mathrm{Broad \, \, Lines}"                 , linestyle=:solid             , linecolor=ith_color(2)                          )

        for (cname, comp) in res.bestfit
            if cname in [:NarrowLines]
                narrow = res.bestfit(:NarrowLines) #Narrow lines components
                plot!(wavecomp , narrow, label=L"\mathrm{Narrow \,\, Lines}", linestyle=:solid, linecolor=ith_color(3))
            end
        end
        
        plot!(wavecomp             , iron              , label=L"\mathrm{Iron}"                              , linestyle=:solid             , linecolor=ith_color(4)                          )

        # plot!(wavecomp             , balmer            , label=L"\mathrm{Balmer}"                              , linestyle=:solid             , linecolor=ith_color(4)                          )
        #plot!(wavecomp             , galax             , label=L"\mathrm{Hos \,\, Galaxy}"                              , linestyle=:solid             , linecolor=ith_color(3)                          )
        plot!(wavecomp             , QSOcont           , label=L"\mathrm{QSO \,\, continuum}"                , linestyle=:solid             , linecolor=ith_color(5)                          )
        plot!(wavecomp             , Nuisance          , label=L"\mathrm{Nuisance \,\, Lines}"               , linestyle=:solid             , linecolor=ith_color(6)                          )
        vline!([656.46]                                , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([486.27]                                , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([279.99]                                , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([1282.0]                                 , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([1083.2]                                 , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([1005.0]                                 , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([500.7]                                  , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([386.9]                                  , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([630.0]                                  , label=""                                            , linestyle=:dash              , color="black"                                   )
        #vline!([671.6]                                  , label=""                                            , linestyle=:dash              , color="black"                                   )
        vline!([232.6]                                  , label=""                                            , linestyle=:dash              , color="black"                                   )
        
        annotate!([650.06]              , [0.3]               , text(L"\mathrm{H\alpha}"                                                    , :black                                      ,12))
        annotate!([615.0]                , [0.3]               , text(L"\mathrm{[OI]\lambda 6300}"                                                         , :black                                      ,12))
        annotate!([681.6]                , [0.05]              , text(L"\mathrm{SII}"                                                        , :black                                      ,12))
        
        if file == -640749119461755401
            annotate!([492.07]          , [4.]                , text(L"\mathrm{H\beta}"                                                     , :black                                      ,12))
            annotate!([515.7]            , [4.]                , text(L"\mathrm{[OIII]\lambda 5007}"                                           , :black                                      ,12))
            annotate!([397.0]            , [4.]                , text(L"\mathrm{[NeIII]}"                                                      , :black                                      ,12))
            
        else
            annotate!([492.07]          , [1.5]               , text(L"\mathrm{H\beta}"                                                     , :black                                      ,12))
        end
        
        annotate!([288.0]                , [30]                , text(L"\mathrm{MgII}"                                                       , :black                                      ,12))
        annotate!([1295.0]               , [0.02]              , text(L"\mathrm{Pa\beta}"                                                    , :black                                      ,12))
        annotate!([1095.0]               , [0.02]              , text(L"\mathrm{HeI}"                                                        , :black                                      ,12))
        annotate!([1020.0]               , [0.02]              , text(L"\mathrm{Pa\delta}"                                                   , :black                                      ,12))
        annotate!([225.0]                , [33]                , text(L"\mathrm{CII] \lambda 2326}"                                                        , :black                                      ,12))
        xlims!(minimum(waveor), maximum(waveor))
        
        if file in [-640749119461755401]
            plot!(legend=false)
            annotate!((0.1,0.9)                               , text(L"\textbf{%$(file)}"                        , :left                    , :black                                      ,12))
                
        elseif file in [-646101801481095171, -603446588505854734] 
            plot!(legend=false)
            annotate!((0.73,0.9)                              , text(L"\textbf{%$(file)}"                        , :left                    , :black                                      ,12))
        else
            plot!(legend=:top)
            annotate!((0.73,0.9)                              , text(L"\textbf{%$(file)}"                        , :left                    , :black                                      ,12))
        end
        
        push!(accum, esp)
        #######
    end

    plot(accum..., layout=grid(length(accum), 1); SERIES_OPTS..., xlabel=LaTeXString("Rest-frame wavelength [nm]"), ylabel=L"$\mathrm{Lum. \, dens. \, [10^{42} \, erg \, s^{-1} \, \AA^{-1}]}$")    
    savefig("$(output)/QSFIT_Ha_Hb_MgII.png")

end

##########################################################################
#Number counts for the paper

@printf("Total number of spectra in Euclid sample: %i \n"                               , nrow(euclid))
@printf("Number of sources in sample with host galaxy contribution: %i \n"              , count(sub.hostgal))
@printf("Number of sources in sample WITHOUT host galaxy contribution: %i \n"           , nrow(euclid)-count(sub.hostgal))
@printf("Number of sources in sample present in Fu et al.: %i \n"                       , count(sub.FU))
@printf("Number of sources in sample present in DESI: %i \n"                            , count(sub.DESI))
@printf("Number of sources in sample present in QUBRICS: %i \n"                         , count(sub.QUBRICS))
@printf("Good quality spectra in Euclid sample: %i \n"                                  , count(qc.good))
@printf("Number of spectra with counterpart in DESI: %i \n"                             , nrow(euclid_match))
@printf("Good quality spectra with counterpart in DESI: %i \n"                          , count(qc_match.good))

@printf("\n")
@printf("--")
@printf("\n")
@printf("Mean H magnitude for quality sample: %.1f +/- %.1f \n", mean(euclid[qc.good, :Hmag]), std(euclid[qc.good, :Hmag]))

@printf("\n")
@printf("--")
@printf("\n")

iii = findall(sub.FU      .& qc.good)
@printf("Good quality sources present in Fu et al.: %i \n", length(iii))
@printf("	redshift range: %.3f -- %.3f \n"          , minimum(euclid[iii, :Redshift]), maximum(euclid[iii, :Redshift]))
@printf("	H mag range: %.2f -- %.2f \n"             , minimum(euclid[iii, :Hmag])    ,  maximum(euclid[iii, :Hmag]))

@printf("\n")

iii = findall(sub.DESI      .& qc.good)
@printf("Good quality sources present in DESI: %i \n", length(iii))
@printf("	redshift range: %.3f -- %.3f \n"     , minimum(euclid[iii, :Redshift]), maximum(euclid[iii, :Redshift]))
@printf("	H mag range: %.2f -- %.2f \n"        , minimum(euclid[iii, :Hmag])    , maximum(euclid[iii, :Hmag]))

@printf("\n")

iii = findall(sub.QUBRICS      .& qc.good)
@printf("Good quality sources present in QUBRICS: %i \n", length(iii))
@printf("	redshift range: %.3f -- %.3f \n"        , minimum(euclid[iii, :Redshift]), maximum(euclid[iii, :Redshift]))
@printf("	H mag range: %.2f -- %.2f \n"           , minimum(euclid[iii, :Hmag])    , maximum(euclid[iii, :Hmag]))

@printf("\n")
@printf("--")
@printf("\n")

@printf("QSO slope for quality sample: %.2f +/- %.2f \n"                  , median(euclid[qc.good, :QSOcont_alpha])                                   , mad(euclid[qc.good, :QSOcont_alpha]                                 , normalize=true))
@printf("Redshift range for quality sample: %.2f -- %.2f \n"              , minimum(euclid[qc.good, :Redshift])                                       , maximum(euclid[qc.good, :Redshift]))
@printf("Bolometric lum. for quality sample: %.2f +/- %.2f \n"            , median(filter(!isnan, log10.(euclid[qc.good, :Lbol_mean])))               , mad(filter(!isnan, log10.(euclid[qc.good, :Lbol_mean]))             , normalize=true))
@printf("\n")
@printf("QSO slope for quality sample at z<0.8: %.2f +/- %.2f \n"         , median(euclid[qc.good .& sub.lowz, :QSOcont_alpha])                       , mad(euclid[qc.good .& sub.lowz, :QSOcont_alpha]                     , normalize=true))
@printf("Redshift range for z<0.8 quality sample: %.2f -- %.2f \n"        , minimum(euclid[qc.good .& sub.lowz, :Redshift])                           , maximum(euclid[qc.good .& sub.lowz, :Redshift]))
@printf("Bolometric lum. for z<0.8 quality sample: %.2f +/- %.2f \n"      , median(filter(!isnan, log10.(euclid[qc.good .& sub.lowz, :Lbol_mean])))   , mad(filter(!isnan, log10.(euclid[qc.good .& sub.lowz, :Lbol_mean])) , normalize=true))
@printf("\n")
@printf("QSO slope for quality sample at 0.8<z<1.9: %.2f +/- %.2f \n"     , median(euclid[qc.good .& sub.cosmo, :QSOcont_alpha])                      , mad(euclid[qc.good .& sub.cosmo, :QSOcont_alpha]                    , normalize=true))
@printf("Redshift range for 0.8<z<1.9 quality sample: %.2f -- %.2f \n"    , minimum(euclid[qc.good .& sub.cosmo, :Redshift])                          , maximum(euclid[qc.good .& sub.cosmo, :Redshift]))
@printf("Bolometric lum. for 0.8<z<1.9 quality sample: %.2f +/- %.2f \n"  , median(filter(!isnan, log10.(euclid[qc.good .& sub.cosmo, :Lbol_mean])))  , mad(filter(!isnan, log10.(euclid[qc.good .& sub.cosmo, :Lbol_mean])), normalize=true))
@printf("\n")
@printf("QSO slope for quality sample z>1.9: %.2f +/- %.2f \n"            , median(euclid[qc.good .& sub.highz, :QSOcont_alpha])                      , mad(euclid[qc.good .& sub.highz, :QSOcont_alpha]                    , normalize=true))
@printf("Redshift range for z>1.9 quality sample: %.2f -- %.2f \n"        , minimum(euclid[qc.good .& sub.highz, :Redshift])                          , maximum(euclid[qc.good .& sub.highz, :Redshift]))
@printf("Bolometric lum. for z>1.9 quality sample: %.2f +/- %.2f \n"      , median(filter(!isnan, log10.(euclid[qc.good .& sub.highz, :Lbol_mean])))  , mad(filter(!isnan, log10.(euclid[qc.good .& sub.highz, :Lbol_mean])), normalize=true))
@printf("\n")
@printf("QSO slope for DESI quality sample: %.2f +/- %.2f \n"             , median(desi[qc_desi.good, :QSOcont_alpha])                                , mad(desi[qc_desi.good, :QSOcont_alpha]                              , normalize=true))
@printf("Redshift range for DESI quality sample: %.2f -- %.2f \n"         , minimum(desi[qc_desi.good, :Redshift])                                    , maximum(desi[qc_desi.good, :Redshift]))
@printf("Bolometric lum. for DESI quality sample: %.2f +/- %.2f \n"       , median(filter(!isnan, log10.(desi[qc_desi.good, :Lbol_mean])))            , mad(filter(!isnan, log10.(desi[qc_desi.good, :Lbol_mean])), normalize=true))

@printf("\n")
@printf("--")
@printf("\n")

@printf("Number of sources with Ha line: %i \n"                      , count(qc.Ha))
@printf("Mean luminosity for Ha: %.2f +/- %.2f \n"                   , mean(log10.(euclid[qc.Ha, :Ha_br_norm]) .+ 42)         , std(log10.(euclid[qc.Ha, :Ha_br_norm]) .+ 42))
@printf("Mean FWHM for Ha: %.2f +/- %.2f \n"                         , mean(log10.(euclid[qc.Ha, :Ha_br_fwhm]))               , std(log10.(euclid[qc.Ha, :Ha_br_fwhm]))) 
@printf("Linear Mean FWHM for Ha: %.2f +/- %.2f \n"                         , mean(euclid[qc.Ha, :Ha_br_fwhm])                , std(euclid[qc.Ha, :Ha_br_fwhm])) 
@printf("Linear Median FWHM for Ha: %.2f +/- %.2f \n"                       , median(euclid[qc.Ha, :Ha_br_fwhm])              , mad(euclid[qc.Ha, :Ha_br_fwhm])) 
@printf("Percentage of sources with Ha FHWM > 10 000 km/s: %.2f %%\n", (count(euclid[qc.Ha, :Ha_br_fwhm] .> 10000.)/length(findall(qc.Ha)))*100)
@printf("Mean v_off for Ha: %.2f +/- %.2f \n"                        , mean(euclid[qc.Ha, :Ha_br_voff])                       , std(euclid[qc.Ha, :Ha_br_voff]))
@printf("Mean BH mass for Ha: %.2f +/- %.2f \n"                      , mean(euclid[qc.Ha, :MBH_Ha_ShenLiu2012])               , std(euclid[qc.Ha, :MBH_Ha_ShenLiu2012]))
@printf("Mean Edd. ratio for Ha (log): %.2f +/- %.2f \n"             , mean(filter(!isnan, log10.(euclid[qc.Ha, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.Ha, :Edd_ratio]))))

@printf("\n")
@printf("--")
@printf("\n")

@printf("Number of sources with Hb line: %i \n"                       , count(qc.Hb))
@printf("Mean luminosity for Hb: %.2f +/- %.2f \n"                    , mean(log10.(euclid[qc.Hb, :Hb_br_norm]) .+ 42)         , std(log10.(euclid[qc.Hb, :Hb_br_norm]) .+ 42))
@printf("Mean FWHM for Hb: %.2f +/- %.2f \n"                          , mean(log10.(euclid[qc.Hb, :Hb_br_fwhm]))               , std(log10.(euclid[qc.Hb, :Hb_br_fwhm]))) 
@printf("Linear Mean FWHM for Hb: %.2f +/- %.2f \n"                         , mean(euclid[qc.Hb, :Hb_br_fwhm])               , std(euclid[qc.Hb, :Hb_br_fwhm]))
@printf("Linear Median FWHM for Hb: %.2f +/- %.2f \n"                       , median(euclid[qc.Hb, :Hb_br_fwhm])              , mad(euclid[qc.Hb, :Hb_br_fwhm])) 
@printf("Percentage of sources with Hb FHWM > 10 000 km/s: %.2f %% \n", (count(euclid[qc.Hb, :Hb_br_fwhm] .> 10000) / count(qc.Hb))*100)

@printf("Mean v_off for Hb: %.2f +/- %.2f \n"           , mean(euclid[qc.Hb, :Hb_br_voff])                       , std(euclid[qc.Hb, :Hb_br_voff]))
@printf("Mean BH mass for Hb: %.2f +/- %.2f \n"         , mean(filter(!isnan, euclid[qc.Hb, :MBH_Hb_WuShen2022])), std(filter(!isnan, euclid[qc.Hb, :MBH_Hb_WuShen2022])))
@printf("Mean Edd. ratio for Hb (log): %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.Hb, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.Hb, :Edd_ratio]))))

@printf("\n")
@printf("--")
@printf("\n")

@printf("Number of sources with MgII line: %i \n"                       , count(qc.MgII))
@printf("Mean luminosity for MgII: %.2f +/- %.2f \n"                    , mean(log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42)  , std(log10.(euclid[qc.MgII, :MgII_2798_br_norm]) .+ 42))
@printf("Mean FWHM for MgII: %.2f +/- %.2f \n"                          , mean(log10.(euclid[qc.MgII, :MgII_2798_br_fwhm]))        ,  std(log10.(euclid[qc.MgII, :MgII_2798_br_fwhm])))
@printf("Percentage of sources with MgII FHWM > 10 000 km/s: %.2f %% \n", (count(euclid[qc.MgII, :MgII_2798_br_fwhm] .> 10000) / count(qc.MgII))*100)

@printf("Mean v_off for MgII: %.2f +/- %.2f \n"           , mean(euclid[qc.MgII, :MgII_2798_br_voff])                , std(euclid[qc.MgII, :MgII_2798_br_voff]))
@printf("Mean BH mass for MgII: %.2f +/- %.2f \n"         , mean(euclid[qc.MgII, :MBH_MgII_WuShen2022])              , std(euclid[qc.MgII, :MBH_MgII_WuShen2022]))
@printf("Mean Edd. ratio for MgII (log): %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.MgII, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.MgII, :Edd_ratio]))))

@printf("\n")
@printf("--")
@printf("\n")

@printf("Number of sources with HeI line: %i \n"                       , count(qc.HeI))
@printf("Mean luminosity for HeI: %.2f +/- %.2f \n"                    , mean(log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42) , std(log10.(euclid[qc.HeI , :HeI_10832_br_norm]) .+ 42))
@printf("Mean FWHM for HeI: %.2f +/- %.2f \n"                          , mean(log10.(euclid[qc.HeI , :HeI_10832_br_fwhm]))       , std(log10.(euclid[qc.HeI , :HeI_10832_br_fwhm])))
@printf("Percentage of sources with HeI FHWM > 10 000 km/s: %.2f %% \n", (count(euclid[qc.HeI, :HeI_10832_br_fwhm] .> 10000) / count(qc.HeI))*100)

@printf("Mean v_off for HeI: %.2f +/- %.2f \n"     , mean(euclid[qc.HeI , :HeI_10832_br_voff])               , std(euclid[qc.HeI , :HeI_10832_br_voff]))
@printf("Mean BH mass for HeI: %.2f +/- %.2f \n"   , mean(euclid[qc.HeI, :MBH_HeI_Ricci])                    , std(euclid[qc.HeI, :MBH_HeI_Ricci]))
@printf("Mean Edd. ratio for HeI (log): %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.HeI, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.HeI, :Edd_ratio]))))

@printf("\n")
@printf("--")
@printf("\n")

@printf("Number of sources with Pab line: %i \n"                       , count(qc.Pab))
@printf("Mean luminosity for Pab: %.2f +/- %.2f \n"                    , mean(log10.(euclid[qc.Pab , :Pab_br_norm]) .+ 42)       , std(log10.(euclid[qc.Pab , :Pab_br_norm]) .+ 42))
@printf("Mean FWHM for Pab: %.2f +/- %.2f \n"                          , mean(log10.(euclid[qc.Pab , :Pab_br_fwhm]))             , std(log10.(euclid[qc.Pab , :Pab_br_fwhm])))
@printf("Percentage of sources with Pab FHWM > 10 000 km/s: %.2f %% \n", (count(euclid[qc.Pab, :Pab_br_fwhm] .> 10000) / count(qc.Pab))*100)

@printf("Mean v_off for Pab: %.2f +/- %.2f \n"           , mean(euclid[qc.Pab , :Pab_br_voff])                     , std(euclid[qc.Pab , :Pab_br_voff]))
@printf("Mean BH mass for Pab: %.2f +/- %.2f \n"         , mean(euclid[qc.Pab, :MBH_Pab_Ricci])                    , std(euclid[qc.Pab, :MBH_Pab_Ricci]))
@printf("Mean Edd. ratio for Pab (log): %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.Pab, :Edd_ratio]))), std(filter(!isnan, log10.(euclid[qc.Pab, :Edd_ratio]))))


@printf("\n")
@printf("--")
@printf("\n")

@printf("Mean bolometric luminosity for good sample: %.2f +/- %.2f \n", mean(filter(!isnan, log10.(euclid[qc.good, :Lbol_mean])))              , std(filter(!isnan, log10.(euclid[qc.good, :Lbol_mean]))))
@printf("Mean BH mass for good sample: %.2f +/- %.2f \n"              , mean(filter(!isnan, euclid[qc.good, :MBH_mean]))                       , std(filter(!isnan, euclid[qc.good, :MBH_mean])))
@printf("Mean Edd. ratio for good sample (log): %.2f +/- %.2f \n"     , mean(filter(!isnan, log10.(euclid[qc.good, :Edd_ratio])))              , std(filter(!isnan, log10.(euclid[qc.good, :Edd_ratio]))))
@printf("Mean Edd. ratio for good sample: %.2f +/- %.2f \n"     , mean(filter(!isnan, euclid[qc.good, :Edd_ratio]))              , std(filter(!isnan, euclid[qc.good, :Edd_ratio])))





