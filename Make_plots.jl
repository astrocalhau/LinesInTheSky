using CSV, Statistics, Plots, StatsPlots, LaTeXStrings, DataFrames, FITSIO, StatsBase, GLM, Printf, Dierckx, QSFit, JSON, SortMerge


#########################################################################################################
#													#
# Script to create the various figures for the Euclid Q1 Paper "Lines in the Sky" 			#
#													#
#########################################################################################################

########################### Make folder for figures #####################################################
mkpath("Figures")

########################### Read CSV table with relevant data ###########################################
############################### and loads it to a DataFrame #############################################
f = FITS("results_Euclid/QSFIT_RESULTS.fits")
df = DataFrame(f[2])
close(f)

#############################################################################
# Identify subsets and quality cut
ii = (FU      = findall(df.FU .== 1),
      DESI    = findall(df.DESI .== 1),
      QUBRICS = findall(df.QUBRICS .== 1),
      hostgal = findall(df.Galaxy_norm .>=0)) #Sources where the host galaxy template was fit (regardless of reliability)

qc = (good = findall( df.good .== 1),                                         #Sources passing the quality cut
      Ha   = findall((df.good .== 1)  .&&  (df.Ha_br_reliable .== 1)),        #quality cut for Ha-based BH masses
      Hb   = findall((df.good .== 1)  .&&  (df.Hb_br_reliable .== 1)),        #quality cut for Hb-based BH masses
      MgII = findall((df.good .== 1)  .&&  (df.MgII_2798_br_reliable .== 1)), #quality cut for MgII-based BH masses
      Pab  = findall((df.good .== 1)  .&&  (df.Pab_br_reliable .== 1)),       #quality cut for Pab-based BH masses
      HeI  = findall((df.good .== 1)  .&&  (df.HeI_10832_br_reliable .== 1))) #quality cut for HeI-based BH masses

palette = [:darkred, :darkgreen, :darkblue, :gray]

#############################################################################
#############################################################################
#############################################################################

# Redshift Histogram
kws = (normalize=false, guidefontsize=16, tickfontsize=16, legendfontsize=14, fillalpha=0.6, alpha=1, grid=true, framestyle=:box, fill=true)
stephist( df[:         , :Redshift], label=L"$\mathrm{Merged \, \, sample}$"                    , bins=30, color=palette[1], fillcolor=palette[1]; kws...)
stephist!(df[ii.FU     , :Redshift], label=L"$\mathrm{Fu \,\, et \,\, al. \,\, (in \,\, prep)}$", bins=30, color=palette[2], fillcolor=palette[2]; kws...)
stephist!(df[ii.DESI   , :Redshift], label=L"$\mathrm{DESI}$"                                   , bins=30, color=palette[3], fillcolor=palette[3]; kws...)
stephist!(df[ii.QUBRICS, :Redshift], label=L"$\mathrm{QUBRICS}$"                                , bins=30, color=palette[4], fillcolor=palette[4]; kws...)
title!(L"$\mathrm{\textbf{Full \,\, sample}}$")
plot!(formatter=:latex)
xlabel!(L"$z$")
ylabel!(L"$\mathrm{Counts}$")
savefig("Figures/Redshift_Hist_Figure.pdf")

#############################################################################
stephist( df[                      qc.good , :Redshift], label=L"$\mathrm{Merged \, \, sample}$"                    , bins=30, color=palette[1], fillcolor=palette[1]; kws...)
stephist!(df[intersect(ii.FU     , qc.good), :Redshift], label=L"$\mathrm{Fu \,\, et \,\, al. \,\, (in \,\, prep)}$", bins=30, color=palette[2], fillcolor=palette[2]; kws...)
stephist!(df[intersect(ii.DESI   , qc.good), :Redshift], label=L"$\mathrm{DESI}$"                                   , bins=30, color=palette[3], fillcolor=palette[3]; kws...)
stephist!(df[intersect(ii.QUBRICS, qc.good), :Redshift], label=L"$\mathrm{QUBRICS}$"                                , bins=30, color=palette[4], fillcolor=palette[4]; kws...)
title!(L"$\mathrm{\textbf{Good\,\, sample}}$")
xlabel!(L"$z$")
ylabel!(L"$\mathrm{Counts}$")
plot!(formatter=:latex)
savefig("Figures/Redshift_Hist_Quality_Figure.pdf")

#############################################################################
#############################################################################
#############################################################################
# All three main BH mass histograms as subplots

function add_Median_Mad(vv)
    mm = median(vv)
    ss = mad(vv, normalize=true)
    smm = @sprintf("%.2f", mm)
    sss = @sprintf("%.2f", ss)
    vline!([median(vv)]               , label=L"$\mathrm{Median} \pm \mathrm{MAD}: %$smm \pm %$sss $", color=:black, linewidth=3)
    vline!( median(vv) .+ [1,-1] .* ss, label=""                                                     , color=:black, linewidth=2, linestyle=:dash)
end

kws = (titlefontsize=20, guidefontsize=20, tickfontsize=16, legendfontsize=16,
       fill=true, color=:royalblue4, grid=true, framestyle=:box, fillalpha=0.8,
       legend=:outertop, background_color_legend=nothing, foreground_color_legend=nothing, legend_column=-1,
       label="", xlabel=L"$\log_{10}(M_{\mathrm{BH}}/\mathrm{M_{\odot}})$", ylabel=L"$\mathrm{Counts}$", xlims=[6,11])

#################################################
# MgII
vv = filter(!isnan, df[qc.MgII, :MBH_MgII_WuShen2022])
stephist(vv, bins=5; kws...)
add_Median_Mad(vv)
annotate!([6.5], [ylims()[2] - (ylims()[2] - ylims()[1]) * 0.1], text(L"$\mathrm{MgII}$", 20, :black , rotation=0))
plot!(formatter=:latex);  savefig("Figures/BH_mass_MgII.pdf")

############################################
# Halpha
vv = filter(!isnan, df[qc.Ha, :MBH_Ha_ShenLiu2012])
stephist(vv, bins=15; kws...)
add_Median_Mad(vv)
annotate!([6.5], [ylims()[2] - (ylims()[2] - ylims()[1]) * 0.1], text(L"$\mathrm{H}\alpha$", 20, :black , rotation=0))
plot!(formatter=:latex);  savefig("Figures/BH_mass_Ha.pdf")

###############################################
# Hbeta
vv = filter(!isnan, df[qc.Hb, :MBH_Hb_WuShen2022])
stephist(vv, bins=10; kws...)
add_Median_Mad(vv)
annotate!([6.5], [ylims()[2] - (ylims()[2] - ylims()[1]) * 0.1], text(L"$\mathrm{H}\beta$", 20, :black , rotation=0))
plot!(formatter=:latex);  savefig("Figures/BH_mass_Hb.pdf")

##################################################
# Pab
vv = filter(!isnan, df[qc.Pab, :MBH_Pab_Ricci])
stephist(vv, bins=5; kws...)
add_Median_Mad(vv)
annotate!([6.5], [ylims()[2] - (ylims()[2] - ylims()[1]) * 0.1], text(L"$\mathrm{Pa}\beta$", 20, :black , rotation=0))
plot!(formatter=:latex);  savefig("Figures/BH_mass_Pab.pdf")

##################################################
# HeI
vv = filter(!isnan, df[qc.HeI, :MBH_HeI_Ricci])
stephist(vv, bins=10; kws...)
add_Median_Mad(vv)
annotate!([6.5], [ylims()[2] - (ylims()[2] - ylims()[1]) * 0.1], text(L"$\mathrm{He\,I}$", 20, :black , rotation=0))
plot!(formatter=:latex);  savefig("Figures/BH_mass_HeI.pdf")

##########################################
# Comparison between Ha and Hb
vv = filter(!isnan, df[intersect(qc.Ha, qc.Hb), :MBH_Ha_ShenLiu2012] .- df[intersect(qc.Ha, qc.Hb), :MBH_Hb_WuShen2022])
stephist(vv, bins=10; kws..., xlims=[-3,3], xlabel=L"$\log_{10}(M_{\mathrm{H}\alpha}/M_{\mathrm{H}\beta})$")
add_Median_Mad(vv)
annotate!([-2], [ylims()[2] - (ylims()[2] - ylims()[1]) * 0.1], text(L"$\mathrm{H}\alpha\,vs\,\mathrm{H}\beta$", 20, :black , rotation=0))
plot!(formatter=:latex);  savefig("Figures/BH_mass_Ha_vs_Hb.pdf")

#############################################################################
#############################################################################
#############################################################################

# QSO Continuum alpha index histogram
#With quality cut and redshift separation
# Creating relevant dataframes for the redshift bins
z1 = df[intersect(qc.good, findall(      df.Redshift .< 1)), :]
z2 = df[intersect(qc.good, findall(1 .<= df.Redshift .< 2)), :]
z3 = df[intersect(qc.good, findall(2 .<= df.Redshift     )), :]

#Estimating medians for the bins
medianz1 = @sprintf("%.2f",median(z1.QSOcont_alpha))
medianz2 = @sprintf("%.2f",median(z2.QSOcont_alpha))
medianz3 = @sprintf("%.2f",median(z3.QSOcont_alpha))
medianzall = @sprintf("%.2f",median(df[qc.good, :QSOcont_alpha]))

#Ploting histograms
kws = (guidefontsize=16, tickfontsize=16, legendfontsize=14, fill=true, fillalpha=0.6, legend=:topright, grid=true, framestyle=:box)
stephist(df[qc.good, :QSOcont_alpha], label=L"$\mathrm{Merged\,sample}$", bins=20, color=palette[4], fillcolor=palette[4]; kws...)
stephist!(z1.QSOcont_alpha          , label=L"$z < 1$"                  , bins=20, color=palette[1], fillcolor=palette[1]; kws...)
stephist!(z2.QSOcont_alpha          , label=L"$1 < z < 2$"              , bins=20, color=palette[2], fillcolor=palette[2]; kws..., z_order=2)
stephist!(z3.QSOcont_alpha          , label=L"$z > 2$"                  , bins=20, color=palette[3], fillcolor=palette[3]; kws...)

#Plotting vertical lines for the medians
# MEDIANIE = vline!([median(filter(!isnan, skipmissing(qualcut.QSOcont_alpha)))], label=L"\mathrm{m}\, \alpha_{\lambda} = %$medianzall", color=:black, linewidth = 2, thickness_scalling =1, linestyle=:solid, z_order=6)
# z11 = vline!([median(filter(!isnan, skipmissing(z1.QSOcont_alpha)))], label=L"\mathrm{m}\, \alpha_{\lambda,\,z<1} = %$medianz1", color=:darkred, linewidth = 3, thickness_scalling =1, linestyle=:dash, z_order=7)
# z22 = vline!([median(filter(!isnan, skipmissing(z2.QSOcont_alpha)))], label=L"\mathrm{m}\, \alpha_{\lambda, \, 1<z<2}= %$medianz2", color=:green, linewidth = 3, thickness_scalling =1, linestyle=:dot, z_order=8)
# z33 = vline!([median(filter(!isnan, skipmissing(z3.QSOcont_alpha)))], label=L"\mathrm{m}\, \alpha_{\lambda,\,z>2} = %$medianz3", color=:midnightblue, linewidth = 3, thickness_scalling =1, linestyle=:dashdot, z_order=9)

plot!(formatter=:latex)
xlabel!(L"$\alpha_{\lambda}$")
ylabel!(L"$\mathrm{Counts}$")
savefig("Figures/QSOcont_alpha_Hist_Redshift_2.pdf")

#############################################################################
#############################################################################
#############################################################################

# Chi2 vs SNR plots
#With Quality cut
scatter(df.redchisq, df.DER_SNR, label=L"$\mathrm{Full \,\, sample}$", mc=palette[1], ms=2, markerstrokewidth=0, dpi=300, guidefontsize=14, tickfontsize=14, legendfontsize=14, grid=false, framestyle=:box)
scatter!(df[qc.good, :redchisq], df[qc.good, :DER_SNR], label=L"$\mathrm{Good \,\, sample}$", mc=:royalblue3, ma=0.6, markershape=:circ, ms=4)
plot!(formatter=:latex, xaxis=:log10, yaxis=:log10)
xlabel!(L"$\mathrm{Reduced} \,\, \chi^2$")
ylabel!(L"$\mathrm{Spectrum\ S/N}$")
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




