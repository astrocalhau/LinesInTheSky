using CSV, Statistics, Plots, StatsPlots, LaTeXStrings, DataFrames, FITSIO, StatsBase, GLM, Printf, Dierckx, QSFit, JSON, SortMerge

######################################################################
# Script to create the various figures for the Euclid Q1 Paper "Lines in the Sky"

mkpath("Figures")

f = FITS("results_Euclid/QSFIT_RESULTS.fits")
df = DataFrame(f[2])
close(f)

# Identify subsets and quality cut
ii = (FU      = findall(df.FU .== 1),
      DESI    = findall(df.DESI .== 1),
      QUBRICS = findall(df.QUBRICS .== 1),
      hostgal = findall(df.Galaxy_norm .>=0)) #Sources where the host galaxy template was fit (regardless of reliability)

qc = (good = findall((df.good .== 1)  .&&  (df.QSOcont_reliable .== 1)),      # sources passing the quality cut
      Ha   = findall((df.good .== 1)  .&&  (df.Ha_br_reliable .== 1)),        # quality cut for Ha-based BH masses
      Hb   = findall((df.good .== 1)  .&&  (df.Hb_br_reliable .== 1)),        # quality cut for Hb-based BH masses
      MgII = findall((df.good .== 1)  .&&  (df.MgII_2798_br_reliable .== 1)), # quality cut for MgII-based BH masses
      Pab  = findall((df.good .== 1)  .&&  (df.Pab_br_reliable .== 1)),       # quality cut for Pab-based BH masses
      HeI  = findall((df.good .== 1)  .&&  (df.HeI_10832_br_reliable .== 1))) # quality cut for HeI-based BH masses

# General plot settings
GEN_OPTS = (formatter=:latex, framestyle=:box, grid=false, thickness_scaling=1.2)
# , titlefontsize=8, guidefontsize=8, tickfontsize=8, legendfontsize=7
# , background_color_legend=nothing, foreground_color_legend=nothing, legend_column=-1
palette = [:darkred, :darkgreen, :darkblue, :gray]
GEN_FILL_ALPHA = 0.6
GEN_ALPHA = 1

xscreen(f) = xlims()[1] + (xlims()[2] - xlims()[1]) * f
yscreen(f) = ylims()[1] + (ylims()[2] - ylims()[1]) * f


#############################################################################
# Redshift histograms

SERIES_OPTS = (bins=0:0.25:maximum(df.Redshift), alpha=GEN_ALPHA, fill=true, fillalpha=GEN_FILL_ALPHA)
plot(; GEN_OPTS..., title=L"$\mathrm{\textbf{Full\ sample}}$", xlabel=L"$z$", ylabel=L"$\mathrm{Counts}$")
stephist!(df[:         , :Redshift], label=L"$\mathrm{Merged\ sample}$"          , color=palette[1]; SERIES_OPTS...)
stephist!(df[ii.FU     , :Redshift], label=L"$\mathrm{Fu\ et\ al.\ (in\ prep.)}$", color=palette[2]; SERIES_OPTS...)
stephist!(df[ii.DESI   , :Redshift], label=L"$\mathrm{DESI}$"                    , color=palette[3]; SERIES_OPTS...)
stephist!(df[ii.QUBRICS, :Redshift], label=L"$\mathrm{QUBRICS}$"                 , color=palette[4]; SERIES_OPTS...)
savefig("Figures/Redshift_Hist_Figure.pdf")

plot(; GEN_OPTS..., title=L"$\mathrm{\textbf{Good\ quality\ sample}}$", xlabel=L"$z$", ylabel=L"$\mathrm{Counts}$")
stephist!(df[                      qc.good , :Redshift], label=L"$\mathrm{Merged\ sample}$"          , color=palette[1]; SERIES_OPTS...)
stephist!(df[intersect(ii.FU     , qc.good), :Redshift], label=L"$\mathrm{Fu\ et\ al.\ (in\ prep.)}$", color=palette[2]; SERIES_OPTS...)
stephist!(df[intersect(ii.DESI   , qc.good), :Redshift], label=L"$\mathrm{DESI}$"                    , color=palette[3]; SERIES_OPTS...)
stephist!(df[intersect(ii.QUBRICS, qc.good), :Redshift], label=L"$\mathrm{QUBRICS}$"                 , color=palette[4]; SERIES_OPTS...)
savefig("Figures/Redshift_Hist_Quality_Figure.pdf")


#############################################################################
# BH mass histograms

function add_details!(vv)
    μ = median(vv)
    σ = mad(vv, normalize=true)
    sμ = @sprintf("%.2f", μ)
    sσ = @sprintf("%.2f", σ)
    vline!([median(vv)]              , label="", color=:black, linewidth=3)
    vline!( median(vv) .+ [1,-1] .* σ, label="", color=:black, linewidth=2, linestyle=:dash)
    annotate!([xscreen(0.03)], [yscreen(0.9)], text(L"$\tilde{\mu}=%$sμ $"   , 12, :black, :left))
    annotate!([xscreen(0.03)], [yscreen(0.8)], text(L"$\tilde{\sigma}=%$sσ $", 12, :black, :left))
end

# Specific plot options
PLOT_OPTS = (legend=:outertop,
             thickness_scaling=1.5, xlabel=L"$\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})$", ylabel=L"$\mathrm{Counts}$")
SERIES_OPTS = (label="", bins=6:0.5:11, color=:royalblue4, alpha=GEN_ALPHA, fill=true, fillalpha=GEN_FILL_ALPHA)

# MgII 
plot(; GEN_OPTS..., PLOT_OPTS..., title=L"$\mathrm{MgII}$")
vv = filter(!isnan, df[qc.MgII, :MBH_MgII_WuShen2022])
stephist!(vv; SERIES_OPTS...)
add_details!(vv)
savefig("Figures/BH_mass_MgII.pdf")

# Halpha
plot(; GEN_OPTS..., PLOT_OPTS..., title=L"$\mathrm{H}\alpha$")
vv = filter(!isnan, df[qc.Ha, :MBH_Ha_ShenLiu2012])
stephist!(vv; SERIES_OPTS...)
add_details!(vv)
savefig("Figures/BH_mass_Ha.pdf")

# Hbeta
plot(; GEN_OPTS..., PLOT_OPTS..., title=L"$\mathrm{H}\beta$")
vv = filter(!isnan, df[qc.Hb, :MBH_Hb_WuShen2022])
stephist!(vv; SERIES_OPTS...)
add_details!(vv)
savefig("Figures/BH_mass_Hb.pdf")

# Pab
plot(; GEN_OPTS..., PLOT_OPTS..., title=L"$\mathrm{Pa}\beta$")
vv = filter(!isnan, df[qc.Pab, :MBH_Pab_Ricci])
stephist!(vv; SERIES_OPTS...)
add_details!(vv)
savefig("Figures/BH_mass_Pab.pdf")

# HeI
plot(; GEN_OPTS..., PLOT_OPTS..., title=L"$\mathrm{He\,I}$")
vv = filter(!isnan, df[qc.HeI, :MBH_HeI_Ricci])
stephist!(vv; SERIES_OPTS...)
add_details!(vv)
savefig("Figures/BH_mass_HeI.pdf")

# Comparison between Ha and Hb
plot(; GEN_OPTS..., PLOT_OPTS..., title=L"$\mathrm{H}\alpha\ vs\ \mathrm{H}\beta$", xlabel=L"$\log_{10}(M_{\mathrm{H}\alpha}/M_{\mathrm{H}\beta})$")
i = intersect(qc.Ha, qc.Hb)
vv = filter(!isnan, df[i, :MBH_Ha_ShenLiu2012] .- df[i, :MBH_Hb_WuShen2022])
stephist!(vv; SERIES_OPTS..., bins=minimum(vv):0.25:maximum(vv))
add_details!(vv)
savefig("Figures/BH_mass_Ha_vs_Hb.pdf")


#############################################################################
# QSO Continuum alpha index histogram

# Redshift bins
i0 = intersect(qc.good)
i1 = intersect(qc.good, findall(      df.Redshift .< 1))
i2 = intersect(qc.good, findall(1 .<= df.Redshift .< 2))
i3 = intersect(qc.good, findall(2 .<= df.Redshift     ))

function add_details!(vv, label, color, SERIES_OPTS; kws...)
    μ = median(vv)
    σ = mad(vv, normalize=true)
    sμ = @sprintf("%.2f", μ)
    sσ = @sprintf("%.2f", σ)
    stephist!(vv, label=label * L"\ (\tilde{\mu}=%$sμ)", color=color; SERIES_OPTS..., kws...)
    vline!([median(vv)], label="", color=color, linewidth=2, linestyle=:dash)
end

PLOT_OPTS = (legend=:topright, xlabel=L"$\alpha_{\lambda}$", ylabel=L"$\mathrm{Counts}$")
SERIES_OPTS = (bins=minimum(df.QSOcont_alpha):0.25:maximum(df.QSOcont_alpha), alpha=GEN_ALPHA, fill=true, fillalpha=GEN_FILL_ALPHA)

#Plot histograms
plot(; GEN_OPTS..., PLOT_OPTS...)
add_details!(df[i0, :QSOcont_alpha], L"$\mathrm{Merged\ sample}$", palette[4], SERIES_OPTS)
add_details!(df[i1, :QSOcont_alpha], L"$z < 1$"                  , palette[1], SERIES_OPTS)
add_details!(df[i2, :QSOcont_alpha], L"$1 < z < 2$"              , palette[2], SERIES_OPTS; z_order=2)
add_details!(df[i3, :QSOcont_alpha], L"$z > 2$"                  , palette[3], SERIES_OPTS)
savefig("Figures/QSOcont_alpha_Hist_Redshift_2.pdf")


#############################################################################
# Chi2 vs SNR plots

PLOT_OPTS = (xlabel=L"$\mathrm{Reduced}\ \chi^2$", ylabel=L"$\mathrm{Spectrum\ S/N}$")

plot(; GEN_OPTS..., PLOT_OPTS...)
scatter!(df[:,       :redchisq], df[:      , :DER_SNR], label=L"$\mathrm{Full\ sample}$"    , mc=palette[1] , ms=2, markerstrokewidth=0)
scatter!(df[qc.good, :redchisq], df[qc.good, :DER_SNR], label=L"$\mathrm{Good \,\, sample}$", mc=:royalblue3, ma=0.6, markershape=:circ, ms=4)
plot!(xaxis=:log10, yaxis=:log10)
savefig("Figures/Chi2vsSNR_cut_Fig.pdf")

#############################################################################
#############################################################################
#############################################################################

# H magnitude w/ source Histogram
#####################################

#Quality cut
mags_cut = @df qualcut stephist(:Hmag, label=L"$\mathrm{Merged\,\,sample}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=50, color=:tomato1, grid=false,framestyle=:box, fill=true, legend=:topleft)
Fumag = @df qualcutFU stephist!(:Hmag, label=L"$\mathrm{Fu \,\, et \,\, al. \,\, (in \,\, prep)}$", normalize=false, guidefontsize=14, tickfontsize=14, legendfontsize=12, alpha=1, bins=30, fill=true, color=:royalblue4, fillcolor=:royalblue4, grid=false, framestyle=:box)
Desmag = @df qualcutDESI stephist!(:Hmag, label=L"$\mathrm{DESI}$", normalize=false, guidefontsize=14, tickfontsize=14, legendfontsize=12, alpha=1, bins=30, fill=true, color=:lightblue1, fillcolor=:lightblue1, grid=false, framestyle=:box)
Qubmag = @df qualcutQUBRICS stephist!(:Hmag, label=L"$\mathrm{QUBRICS}$", normalize=false, guidefontsize=16, tickfontsize=16, legendfontsize=12, alpha=1, bins=30, fill=true, color=:skyblue2, fillcolor=:skyblue2, grid=false, framestyle=:box)
title!(L"$\mathrm{\textbf{Good\,\, sample}}$")

plot!(formatter=:latex)
xlabel!(L"$H_E \, \, \mathrm{magnitude}$")
ylabel!(L"$\mathrm{Counts}$")

savefig(mags_cut, "Figures/Hmag_Hist_source_Fig.pdf")

######################################

#Hmag vs Redshift With Quality cut

HmagqualFU = @df qualcutFU scatter(:Redshift, :Hmag, label=L"$\mathrm{Fu \, \, et \, \, al.}$", mc=:royalblue4, markershape=:rect, ms=6, ma=1, dpi=300, guidefontsize=14, tickfontsize=14, legendfontsize=14, grid=false, framestyle=:box)
HmagqualDESI = @df qualcutDESI scatter!(:Redshift, :Hmag, label=L"$\mathrm{DESI}$", mc=:lightblue1, markershape=:pentagon, ms=6, ma=1, dpi=300, guidefontsize=14, tickfontsize=14, legendfontsize=14, grid=false, framestyle=:box)
HmagqualQUBRICS = @df qualcutQUBRICS scatter!(:Redshift, :Hmag, label=L"$\mathrm{QUBRICS}$", mc=:skyblue2, markershape=:utriangle, ms=6, ma=1, dpi=300, guidefontsize=16, tickfontsize=16, legendfontsize=14, grid=false, framestyle=:box, legend=:bottomright)

title!(L"$\mathrm{\textbf{Good \,\, sample}}$")
ylims!(14,23)
plot!(formatter=:latex)
xlabel!(L"$z$")
ylabel!(L"$H_E$")

savefig(HmagqualFU, "Figures/Hmag_vs_z_cut_Fig.pdf")


#############################################################################
#############################################################################
#############################################################################

# Luminosity Histograms
########################################
# Ha Histogram

dfgoodHa=df[coalesce.(df.Ha_br_reliable .> 0, false), :]
totalha = log10.((dfgoodHa.Ha_br_norm).*1E42)
dfgoodHacut=qualcut[coalesce.(qualcut.Ha_br_reliable .> 0, false), :]
totalhacut = log10.((dfgoodHacut.Ha_br_norm).*1E42)


#Quality cut
Ha_cut = @df dfgoodHacut stephist(totalhacut, label="", guidefontsize=18, tickfontsize=18, legendfontsize=16, bins=30, fill=true, color=:royalblue3, legend=:topleft, grid=false, framestyle=:box, bottom_margin=2*Plots.mm)

#title!(L"$\mathrm{\textbf{Good\,\, sample}}$")
plot!(formatter=:latex)
xlabel!(L"$\log_{10}(L_{\mathrm{H\alpha}}/\mathrm{erg \, s^{-1}})$")
ylabel!(L"$\mathrm{Counts}$")

savefig(Ha_cut, "Figures/Ha_Lum_Fig.pdf")

#########################################################################################################
#Pab, HeI, Hb Luminosity Histograms
#Creating the quantities for the different lines
GoodHb = df[coalesce.(df.Hb_br_reliable .>0, false), :] #all sources passing the reliability cut
GoodHbcut = qualcut[coalesce.(qualcut.Hb_br_reliable .>0), :] #sources passing the quality cut and reliability cut
totalhb = log10.((GoodHb.Hb_br_norm).*1E42)
totalhbcut = log10.((GoodHbcut.Hb_br_norm).*1E42)

GoodPab = df[coalesce.(df.Pab_br_reliable .>0, false), :]
GoodPabcut = qualcut[coalesce.(qualcut.Pab_br_reliable .>0), :]
totalpab = log10.((GoodPab.Pab_br_norm).*1E42)
totalpabcut = log10.((GoodPabcut.Pab_br_norm).*1E42)

GoodHeI = df[coalesce.(df.HeI_10832_br_reliable .>0, false), :]
GoodHeIcut = qualcut[coalesce.(qualcut.HeI_10832_br_reliable .>0), :]
totalhei = log10.((GoodHb.HeI_10832_br_norm).*1E42)
totalheicut = log10.((GoodHeIcut.HeI_10832_br_norm).*1E42)

GoodMgII = df[coalesce.(df.MgII_2798_br_reliable .>0, false), :]
GoodMgIIcut = qualcut[coalesce.(qualcut.MgII_2798_br_reliable .>0), :]
totalmgii = log10.((GoodMgII.MgII_2798_br_norm).*1E42)
totalmgiicut = log10.((GoodMgIIcut.MgII_2798_br_norm).*1E42)


#Quality cut
Hb_cut = @df GoodHbcut stephist(totalhbcut, label=L"$\mathrm{H\beta}$", guidefontsize=17, tickfontsize=17, legendfontsize=15, bins=10, fill=true, color=:royalblue4, legend=:topleft, grid=false, framestyle=:box, bottom_margin=2*Plots.mm, normalize=:probability)

HeI_cut = @df GoodHeIcut stephist!(totalheicut, label=L"$\mathrm{He\,I}$", guidefontsize=17, tickfontsize=17, legendfontsize=15, bins=5, fill=true, color=:skyblue2, legend=:topleft, grid=false, framestyle=:box, alpha=0.8, bottom_margin=2*Plots.mm, normalize=:probability)

Pab_cut = @df GoodPabcut stephist!(totalpabcut, label=L"$\mathrm{Pa}\beta$", guidefontsize=17, tickfontsize=17, legendfontsize=15, bins=15, fill=true, color=:lightblue1, legend=:topleft, grid=false, framestyle=:box, bottom_margin=2*Plots.mm, normalize=:probability, alpha=0.8)

MgII_cut = @df GoodMgIIcut stephist!(totalmgiicut, label=L"$\mathrm{Mg\,II}$", guidefontsize=17, tickfontsize=17, legendfontsize=15, bins=5, fill=true, color=:tomato1, legend=:topleft, grid=false, framestyle=:box, bottom_margin=2*Plots.mm, normalize=:probability, alpha=0.7)

plot!(formatter=:latex)
xlabel!(L"$\log_{10}(L_{\mathrm{line}}/\mathrm{erg \, s^{-1}})$")

savefig(Hb_cut, "Figures/Hb_Lum_Fig.pdf")

#############################################################################
#############################################################################
#############################################################################

#Velocity offsets and FWHM
#relevant cuts

GoodMgII = df[coalesce.(df.MgII_2798_br_reliable .> 0, false), :]
GoodMgIIcut = qualcut[coalesce.(qualcut.MgII_2798_br_reliable .> 0), :]

badvoffHa = dfgoodHacut[coalesce.(dfgoodHacut.Ha_br_voff .> 450),:]
badvoffHb = GoodHbcut[coalesce.(GoodHbcut.Hb_br_voff .> 450),:]
badvoffMgII = GoodMgIIcut[coalesce.(GoodMgIIcut.MgII_2798_br_voff .> 450),:]
badvoffPab = GoodPabcut[coalesce.(GoodPabcut.Pab_br_voff .> 450),:]
badvoffHeI = GoodHeIcut[coalesce.(GoodHeIcut.HeI_10832_br_voff .> 450),:]

#FWHM Histogram panel
#################################################################

# Ha FWHM Histogram


#Quality cut
fwhm_cut = @df dfgoodHacut stephist(:Ha_br_fwhm, label=L"$\mathrm{\textbf{H\alpha}}$", guidefontsize=18, tickfontsize=18, legendfontsize=18, bins=20, fill=true, color=:royalblue3, grid=false,framestyle=:box)#, title=L"$\mathrm{\textbf{H\alpha}}$")

fwhm_badvoff = @df badvoffHa stephist!(:Ha_br_fwhm, label=L"$\mathrm{\textbf{v_{off}>450\,km\,s^{-1}}}$", guidefontsize=18, tickfontsize=18, legendfontsize=18, bins=20, fill=true, color=:tomato1, grid=false,framestyle=:box)

#title!("# of datapoints used by QSFIT")
plot!(formatter=:latex)
xlabel!(L"$ \mathrm{FWHM \, (km \, s^{-1})}$")
ylabel!(L"$\mathrm{Counts}$")
stephist!(xticks!([0, 5000, 10000, 15000],[L"0", L"5000", L"10\,000", L"15\,000"]))


#################################################################
#################################################################
# Hb FWHM Histogram


#Quality cut
fwhm_cutHb = @df GoodHbcut stephist(:Hb_br_fwhm, label=L"$\mathrm{\textbf{H\beta}}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=20, fill=true, color=:royalblue3, grid=false,framestyle=:box, legend=:topright)#, title=L"$\mathrm{\textbf{H\beta}}$")

fwhm_badvoffHb = @df badvoffHb stephist!(:Hb_br_fwhm, label=L"$\mathrm{\textbf{v_{off}>450\,km\,s^{-1}}}$", guidefontsize=18, tickfontsize=18, legendfontsize=18, bins=20, fill=true, color=:tomato1, grid=false,framestyle=:box)

#title!("# of datapoints used by QSFIT")
plot!(formatter=:latex)
xlabel!(L"$ \mathrm{FWHM \, (km \, s^{-1})}$")
ylabel!(L"$\mathrm{Counts}$")
stephist!(xticks!([0, 5000, 10000, 15000],[L"0", L"5000", L"10\,000", L"15\,000"]))
xlims!(0,15000)

#################################################################
#################################################################
# Mg II FWHM Histogram

#Quality cut
fwhm_cutMg = @df GoodMgIIcut stephist(:MgII_2798_br_fwhm, label=L"$\mathrm{\textbf{Mg\,II}}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=5, fill=true, color=:royalblue3, grid=false,framestyle=:box, legend=:topright)#,  title=L"$\mathrm{\textbf{Mg\,II}}$")

fwhm_badvoffMgII = @df badvoffMgII stephist!(:MgII_2798_br_fwhm, label=L"$\mathrm{\textbf{v_{off}>450\,km\,s^{-1}}}$", guidefontsize=18, tickfontsize=18, legendfontsize=18, bins=5, fill=true, color=:tomato1, grid=false,framestyle=:box)

stephist!(xticks!([0, 5000, 10000, 15000],[L"0", L"5000", L"10\,000", L"15\,000"]))
xlims!(0,15050)
plot!(formatter=:latex)
xlabel!(L"$ \mathrm{FWHM \, (km \, s^{-1})}$")
ylabel!(L"$\mathrm{Counts}$")

################################################################
################################################################
# Pab FWHM Histogram

#Quality cut
fwhm_cutPab = @df GoodPabcut stephist(:Pab_br_fwhm, label=L"$\mathrm{\textbf{Pa\beta}}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=10, fill=true, color=:royalblue3, grid=false,framestyle=:box, legend=:topright)#,  title=L"$\mathrm{\textbf{Mg\,II}}$")

fwhm_badvoffPab = @df badvoffPab stephist!(:Pab_br_fwhm, label=L"$\mathrm{\textbf{v_{off}>450\,km\,s^{-1}}}$", guidefontsize=18, tickfontsize=18, legendfontsize=18, bins=20, fill=true, color=:tomato1, grid=false,framestyle=:box)

stephist!(xticks!([0, 5000, 10000, 15000],[L"0", L"5000", L"10\,000", L"15\,000"]))
xlims!(0,15050)
plot!(formatter=:latex)
xlabel!(L"$ \mathrm{FWHM \, (km \, s^{-1})}$")
ylabel!(L"$\mathrm{Counts}$")

################################################################
################################################################
# Mg II FWHM Histogram

#Quality cut
fwhm_cutHeI = @df GoodHeIcut stephist(:HeI_10832_br_fwhm, label=L"$\mathrm{\textbf{He\,I}}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=10, fill=true, color=:royalblue3, grid=false,framestyle=:box, legend=:topright)#,  title=L"$\mathrm{\textbf{Mg\,II}}$")

fwhm_badvoffHeI = @df badvoffHeI stephist!(:HeI_10832_br_fwhm, label=L"$\mathrm{\textbf{v_{off}>450\,km\,s^{-1}}}$", guidefontsize=18, tickfontsize=18, legendfontsize=18, bins=20, fill=true, color=:tomato1, grid=false,framestyle=:box)

stephist!(xticks!([0, 5000, 10000, 15000],[L"0", L"5000", L"10\,000", L"15\,000"]))
xlims!(0,15050)
plot!(formatter=:latex)
xlabel!(L"$ \mathrm{FWHM \, (km \, s^{-1})}$")
ylabel!(L"$\mathrm{Counts}$")

fwhm_panel = plot(fwhm_cutHb, fwhm_cutMg, fwhm_cut, fwhm_cutPab, fwhm_cutHeI, layout=grid(2, 3, widths=(1/3, 1/3, 1/3)), size=(3600, 1800), margin=15*Plots.mm,right_margin=22*Plots.mm, left_margin=25*Plots.mm, bottom_margin=30*Plots.mm, titlefontsize=47,guidefontsize=45, tickfontsize=45, legendfontsize=35)

savefig(fwhm_panel, "Figures/FWHM_Panel.pdf")

#############################################################################
#############################################################################
#############################################################################
# Bol Lum Mean Histogram
#Defining bins
cosmorange= df[df.Redshift .> 0.9 .&& df.Redshift .< 1.8, :]
cosmorangequal = qualcut[qualcut.Redshift .> 0.9 .&& qualcut.Redshift .< 1.8, :]
qualcut = df[(df.good .== 1), :]

#Quality cut
#masses histogram
massescutMean = @df qualcut stephist(:MBH_mean, label=L"$\mathrm{Good\,\,sample}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=30, fill=true, color=:royalblue3, grid=false, framestyle=:box, xlabel=L"$\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})$", ylabel=L"$\mathrm{Counts}$")

Cosmomean = @df cosmorangequal stephist!(:MBH_mean, label=L"$0.8<z<1.9$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=30, fill=true, color=:lightblue2, fillcolor=:lightblue2, grid=false, framestyle=:box, legend=:topleft)

#statistics
MEANFull = mean(filter(!isnan, skipmissing(qualcut.MBH_mean)))
Meanline = vline!([MEANFull], label=L"$\mathrm{Geo. \, \,Mean}$", color="red", linewidth = 3.5, thickness_scalling =1, linestyle=:dash)


#Luminosities histogram
bolcutMean = @df qualcut stephist(log10.(:Lbol_mean), label=L"$\mathrm{Good\,\,sample}$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=30, fill=true, color=:royalblue3, grid=false, framestyle=:box, xlabel=L"$\log_{10}(L_{bol}/\mathrm{erg \, s^{-1}})$", ylabel=L"$\mathrm{Counts}$")

Cosmomeanlum = @df cosmorangequal stephist!(log10.(:Lbol_mean), label=L"$0.8<z<1.9$", guidefontsize=14, tickfontsize=14, legendfontsize=14, bins=10, fill=true, color=:lightblue2, fillcolor=:lightblue2, grid=false, framestyle=:box)

		#statistics
MEANFull = mean(filter(!isnan, skipmissing(log10.(qualcut.Lbol_mean))))
Meanline = vline!([MEANFull], label=L"$\mathrm{Geo. \, \,Mean}$", color="red", linewidth = 3.5, thickness_scalling =1, linestyle=:dash)
MEDIANFull = median(filter(!isnan, skipmissing(qualcut.MBH_mean)))
MADFull = mad(filter(!isnan, skipmissing(qualcut.MBH_mean)), normalize=true)

plot!(formatter=:latex)



Lummasspanel = plot(bolcutMean, massescutMean, layout=grid(1, 2, widths=(4/8, 4/8)), size=(1600,600), margin=5*Plots.mm, left_margin=10*Plots.mm, bottom_margin=12*Plots.mm , guidefontsize=22, tickfontsize=22, legendfontsize=20, titlefontsize=22)

savefig(Lummasspanel, "Figures/Mean_Bol_Mass_fig.pdf")



#############################################################################
#############################################################################
#############################################################################

######################### WITH DESI #########################################

#############################################################################
#############################################################################
#############################################################################

#Read the DESI results catalog from "run_DESI.jl"
DataDESI = FITS("results_DESI/QSFIT_RESULTS.fits")
dfDESI = DataFrame(DataDESI[2])
close(DataDESI)

######################## MATCHING DESI catalogue WITH EUCLID catalogue ######
j = sortmerge(df.ID, dfDESI.ID_EUCLID)
dfmatched = df[j[1],:]
dfDESImatched = dfDESI[j[2],:]

# Restrictions on the FWHM
dffwhm = dfmatched[coalesce.(dfmatched.Ha_br_fwhm .> 2000., false) .&& coalesce.(dfmatched.Ha_br_fwhm .< 15000.,false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .> 2000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .< 15000.,false),:]

dffwhm_qual = dfmatched[coalesce.(dfmatched.Ha_br_fwhm .> 2000., false) .&& coalesce.(dfmatched.Ha_br_fwhm .< 15000.,false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .> 2000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .< 15000.,false) .&& coalesce.(dfmatched.good .== 1., false) .&& coalesce.(dfDESImatched.good .== 1., false),:]


#Quality cut with FWHM restrictions
# For Euclid data
qualcut = dfmatched[coalesce.(dfmatched.good .== 1, false) .&& coalesce.(dfDESImatched.good .== 1, false) .&& coalesce.(dfmatched.Ha_br_fwhm .> 2000., false) .&& coalesce.(dfmatched.Ha_br_fwhm .< 15000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .> 2000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .< 15000.,false),:]

#For DESI data
qualcutDESI = dfDESImatched[coalesce.(dfmatched.good .== 1, false) .&& coalesce.(dfDESImatched.good .== 1, false) .&& coalesce.(dfmatched.Ha_br_fwhm .> 2000., false) .&& coalesce.(dfmatched.Ha_br_fwhm .< 15000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .> 2000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .< 15000.,false),:]


HaFWHM_cut = dfmatched[(dfmatched.good .== 1) .&& coalesce.(dfmatched.Ha_br_fwhm .>2000, false) .&& coalesce.(dfmatched.Ha_br_fwhm .<15000, false),:] #coalesce is needed for


######################################################################################
## DESI BH Mass - Euclid Ha BH mass
minusHa = dfmatched.MBH_Ha_ShenLiu2012 .-dfDESImatched.MBH_MgII_WuShen2022 #difference between Ha Euclid BHM and MgII DESI BHM
minusHacut = qualcut.MBH_Ha_ShenLiu2012 .-qualcutDESI.MBH_MgII_WuShen2022 #same but with qualitry cut applied
minusHaFWHM = qualcut.MBH_Ha_ShenLiu2012 .-qualcutDESI.MBH_MgII_WuShen2022



MBHDESIHaEuclidcut = @df HaFWHM_cut stephist(minusHaFWHM, xlabel=L"$\log_{10}(M_{\mathrm{BH},\, Euclid}/\mathrm{M_{\odot}})-\log_{10}(M_{\mathrm{BH},\, \mathrm{DESI}}/\mathrm{M_{\odot}})$", label=L"$\mathrm{H}\alpha - \mathrm{MgII}$", guidefontsize=12, tickfontsize=12, legendfontsize=12, bins=10, fill=true, color=:royalblue3, grid=false, legend=:topleft, framestyle=:box)

madnorm = mad(filter(!isnan, skipmissing(minusHaFWHM)), normalize=true)
MEDIAN = median(filter(!isnan, skipmissing(minusHaFWHM)))

number = @sprintf("%.2f", MEDIAN)
madnumber = @sprintf("%.2f", madnorm)

mmean = mean(filter(!isnan, skipmissing(minusHaFWHM)))

#meanline = vline!([mmean], label=L"$\mathrm{mean}$", color="red", linewidth = 3, thickness_scalling =1, linestyle=:dash, z_order=5)
annotate!([0.65], [30], text(L"$\mathrm{Median}=%$number $",20, :black , rotation=0))
annotate!([0.65], [27], text(L"$\mathrm{MAD}=%$madnumber $",20, :black , rotation=0))

plot!(formatter=:latex)

#################################################################
#################################################################
betterqualcut = dfmatched[coalesce.(dfmatched.good .== 1, false) .&& coalesce.(dfDESImatched.good .== 1, false) .&& coalesce.(dfmatched.Hb_br_fwhm .> 2000., false) .&& coalesce.(dfmatched.Hb_br_fwhm .< 15000., false) .&& coalesce.(dfmatched.Hb_br_norm .> 0., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .> 2000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .< 15000.,false) .&& coalesce.(dfDESImatched.MgII_2798_br_norm .> 0.,false),:]

betterqualcutDESI = dfDESImatched[coalesce.(dfmatched.good .== 1, false) .&& coalesce.(dfDESImatched.good .== 1, false) .&& coalesce.(dfmatched.Hb_br_fwhm .> 2000., false) .&& coalesce.(dfmatched.Hb_br_fwhm .< 15000., false) .&& coalesce.(dfmatched.Hb_br_norm .> 0., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .> 2000., false) .&& coalesce.(dfDESImatched.MgII_2798_br_fwhm .< 15000.,false) .&& coalesce.(dfDESImatched.MgII_2798_br_norm .> 0.,false),:]

MgIIoverHb = betterqualcutDESI.MgII_2798_br_norm ./ betterqualcut.Hb_br_norm

MgIIHbhist = stephist(MgIIoverHb, xlabel=L"$\mathrm{L}_{\mathrm{Mg\,II},\, \mathrm{DESI}}(\mathrm{erg \, s^{-1}})/\mathrm{L}_{\mathrm{H\beta},\, EUCLID} (\mathrm{erg \, s^{-1}})$", label=L"$\mathrm{H}\beta / \mathrm{MgII}$", guidefontsize=12, tickfontsize=12, legendfontsize=12, bins=20, fill=true, color=:white, fillcolor=:royalblue1, grid=false, legend=:topright, framestyle=:box)

madnorm = mad(filter(!isnan, skipmissing(MgIIoverHb)), normalize=true)
MEDIAN = median(filter(!isnan, skipmissing(MgIIoverHb)))

numbermedian = @sprintf("%.1f", MEDIAN)
madnumber = @sprintf("%.1f", madnorm)

mmean = mean(filter(!isnan, skipmissing(MgIIoverHb)))
ratiomean = @sprintf("%.1f", mmean)


Meanline = vline!([mmean], label=L"$\mathrm{mean}= %$ratiomean $", color="red", linewidth = 3.5, thickness_scalling =1, linestyle=:dash)
Medianline = vline!([MEDIAN], label=L"$\mathrm{median}= %$numbermedian $", color=:chartreuse2, linewidth = 3.5, thickness_scalling =1, linestyle=:dot)


plot!(formatter=:latex)
savefig(MgIIHbhist, "Figures/MgIIDESI_HbEuclid_ratio.pdf")

#################################################################
# BHM scatter
jj = sortmerge(qualcut.ID, qualcutDESI.ID_EUCLID)
BHEumatch = qualcut[jj[1],:]
BHDESmatch = qualcutDESI[jj[2],:]

#Quality cut
BHM_cutMg = scatter(BHDESmatch.MBH_MgII_WuShen2022, BHEumatch.MBH_Ha_ShenLiu2012, label="", guidefontsize=14, tickfontsize=14, markershape=:square, legendfontsize=14, bins=50, color=:royalblue3, grid=false,framestyle=:box, legend=:topleft)

identity2 = @df df plot!(:MBH_Ha_ShenLiu2012, :MBH_Ha_ShenLiu2012, label = L"$1:1$", linecolor=:black, ls=:dash, lw=2, z_order=1)

qualcut.MBH_Ha_ShenLiu2012 = replace(qualcut.MBH_Ha_ShenLiu2012, NaN=>missing)
qualcutDESI.MBH_MgII_WuShen2022 = replace(qualcutDESI.MBH_MgII_WuShen2022, NaN=>missing)

plot!(formatter=:latex)
xlabel!(L"$\log_{10}[M_{\mathrm{BH},\,\mathrm{DESI}}(MgII)/\mathrm{M_{\odot}}]$")
ylabel!(L"$\log_{10}[M_{\mathrm{BH},\, Euclid}(\mathrm{H}\alpha)/\mathrm{M_{\odot}}]$")

xlims!(7.5,10)
ylims!(7.5,10)

#savefig(BHM_cutMg, "Figures/BHM_Euclid_DESI_scatter_Fig.png")


panelfwhm = plot(MBHDESIHaEuclidcut, BHM_cutMg ,layout=grid(1, 2, widths=(4/8, 4/8)), size=(1600,600), margin=5*Plots.mm, left_margin=20*Plots.mm, bottom_margin=15*Plots.mm, titlefontsize=23,guidefontsize=23, tickfontsize=23, legendfontsize=23)
savefig(panelfwhm, "Figures/BHM_Euclid_DESI_Panel_Fig.pdf")

##################################################################################
# Chi2/SNR scatter plot
qualcut_DESI = dfDESI[(dfDESI.good .== 1) .&& (dfDESI.redchisq .<6), :]

#With Quality cut
chicut_log = log10.(qualcut_DESI.redchisq)
Kamehameha = @df dfDESI scatter(:redchisq, :DER_SNR, label=L"$\mathrm{Full \,\, sample}$", mc=:salmon, ms=2, markerstrokewidth=0, ma=1, dpi=300, guidefontsize=14, tickfontsize=14, legendfontsize=14, grid=false, framestyle=:box)
Kamecut = @df qualcut_DESI scatter!(:redchisq, :DER_SNR, label=L"$\mathrm{Good\,\,sample}$", mc=:royalblue3, markershape=:rect, ms=4, ma=1, dpi=300, guidefontsize=14, tickfontsize=14, legendfontsize=14, grid=false, framestyle=:box, legend=:bottomright)

title!(L"$\mathrm{\textbf{DESI \,\, spectra}}$")
plot!(formatter=:latex, xaxis=:log10, yaxis=:log10)
plot!(xticks=([1, 2, 5, 10, 20, 50, 100, 500, 2000, 10000],[L"$1$",L"$2$",L"$5$", L"$10$",L"$20$",L"$50$", L"$100$", L"$500$", L"$2000$", L"$10000$"]), yticks=([2, 5, 10, 20],[L"$2$",L"$5$", L"$10$",L"$20$"]))
xlabel!(L"$\mathrm{reduced}\,\, \chi^2$")
ylabel!(L"$\mathrm{S/N_{spectrum}}$")

savefig(Kamehameha, "Figures/Chi2vsSNR_cut_Fig_DESI.pdf")

######################################################################################

#############################################################################
#############################################################################
#############################################################################

# QSO Continuum alpha index histogram
DESI_qualcut = dfDESImatched[coalesce.(dfDESImatched.good .==1, false),:] #for the QSO continuum slope histogram
Dalpha_qual = DESI_qualcut.QSOcont_alpha
qualcut = df[(df.good .== 1),:] #quality cut for the Euclid data


Dalphacut = @df DESI_qualcut stephist(:QSOcont_alpha, label=L"$\mathrm{DESI \,\,sample}$", bins=20, guidefontsize=16, tickfontsize=16, legendfontsize=14, fill=true, color=:royalblue3, legend=:topright, grid=false, framestyle=:box, z_order=1, normalize=:probability)

alphacut = @df qualcut stephist!(:QSOcont_alpha, label=L"$Euclid\,\, \mathrm{sample}$", bins=80, guidefontsize=16, tickfontsize=16, legendfontsize=14, fill=false, color=:black, fillcolor=:black, legend=:topright, grid=false, framestyle=:box,z_order=2, normalize=:probability)

#Estimating the median for the DESI slope
medianDESI  = @sprintf("%.2f",median(filter(!isnan, skipmissing(Dalpha_qual))))

#Drawing the vertical lines for the medians
z15 = vline!([median(filter(!isnan, skipmissing(Dalpha_qual)))], label=L"$\mathrm{m}\,\alpha_{\lambda,\,\mathrm{DESI}}(0.35<z<2.5)=-1.66$", color=:lightsalmon, linewidth = 5, thickness_scalling =1, linestyle=:dash, z_order=3,legendfontsize=14)

z33 = vline!([median(filter(!isnan, skipmissing(z3.QSOcont_alpha)))], label=L"$\mathrm{m}\,\alpha_{\lambda,\,Euclid}(z>2)=-1.71$", color=:midnightblue, linewidth = 3, thickness_scalling =1, linestyle=:dashdot, z_order=4,legendfontsize=12)

plot!(formatter=:latex)
xlabel!(L"$\alpha_{\lambda}$")
xlims!(-5,6)

savefig(Dalphacut, "Figures/QSOcont_alpha_DESI_Fig.pdf")

#############################################################################

#################################################################
#################################################################
#################################################################
#			SDSS plots 				#
#################################################################
#################################################################
#################################################################

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


MBHMgIIZcut = @df qualcutMgII scatter!(:Redshift, :MBH_MgII_WuShen2022, label=L"$\mathrm{Mg \, II}$", mc=:black, markershape=:utriangle, ms=6, ma=1, dpi=300, guidefontsize=18, tickfontsize=18, legendfontsize=18, grid=false, legend=:bottomright,framestyle=:box)

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
xlabel!(L"$\log_{10}(L_{\mathrm{bol}}/\mathrm{erg \, s^{-1}})$")
ylabel!(L"$\log_{10}(M_{\mathrm{BH}}/{\mathrm{M_{\odot}}})$")
ylims!(4,12)
xlims!(42.5,49)

#savefig(SDSSMgIIbol, "MBHvLbol_w_SDSS.png")
#savefig(SDSSMgIIbol, "MBHvLbol_w_SDSS.pdf")
AllSDSS = plot(MBHPabZcut, MBHHaZcutbol, layout=grid(1, 2, widths=(4/8, 4/8)), size=(1600,600), margin=5*Plots.mm, left_margin=10*Plots.mm, bottom_margin=13*Plots.mm, titlefontsize=22, guidefontsize=22, tickfontsize=22, legendfontsize=20, annotationfontsize=20)
plot!(formatter=:latex)
savefig(AllSDSS, "Figures/MBH_SDSS_all.png")
