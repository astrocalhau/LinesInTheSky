using Revise, FITSIO, DataFrames, Gnuplot, Dierckx, ProgressMeter, Statistics, StatsBase, Serialization, CSV
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer, CMPFit

Gnuplot.options.term = "qt size 1600,900 enhanced font 'Latin Modern Roman, 13' lw 1.5"

# ====================================================================
struct SingleSpec
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
    z::Float64
end


abstract type AbstractCompositeBin end

struct LinCompositeBin <: AbstractCompositeBin
    wavelength::Float64
    ispec::Vector{Int64}
    ref::Vector{Float64}
    scaled::Vector{Float64}
end
LinCompositeBin(l::Float64) = LinCompositeBin(l, Vector{Int64}(), Vector{Float64}(), Vector{Float64}())

struct LogCompositeBin <: AbstractCompositeBin
    wavelength::Float64
    ispec::Vector{Int64}
    ref::Vector{Float64}
    scaled::Vector{Float64}
end
LogCompositeBin(l::Float64) = LogCompositeBin(l, Vector{Int64}(), Vector{Float64}(), Vector{Float64}())

LinCompositeBin(cc::LogCompositeBin) = LinCompositeBin(cc.wavelength, cc.ispec,  10 .^(cc.ref),  10 .^(cc.scaled))
LogCompositeBin(cc::LinCompositeBin) = LogCompositeBin(cc.wavelength, cc.ispec, log10.(cc.ref), log10.(cc.scaled))


import Statistics: mean, std
mean(bin::AbstractCompositeBin) =   mean(bin.scaled)
std( bin::AbstractCompositeBin) =   std( bin.scaled)
nn(  bin::AbstractCompositeBin) = length(bin.scaled)

function add_spec!(bin::LinCompositeBin, ispec::Int64, ref::Float64)
    @assert !isnan(ref)
    push!(bin.ispec, ispec)
    push!(bin.ref, ref)
    push!(bin.scaled, ref)
end

save_scaled_as_ref!(bin::AbstractCompositeBin) = bin.ref .= bin.scaled

apply_scale!(bin::LinCompositeBin, scales::Vector{Float64}) = bin.scaled .= bin.ref .* scales[bin.ispec]

function apply_scale!(bin::LinCompositeBin, scale::Float64, ispec::Int)
    i = findfirst(bin.ispec .== ispec)
    bin.scaled[i] = bin.ref[i] * scale
end

function getscaled(bin::AbstractCompositeBin, ispec::Int)
    i = findfirst(bin.ispec .== ispec)
    isnothing(i)  &&  (return (nothing, nothing))
    return (bin.wavelength, bin.scaled[i])
end

function getscaled(bins::Vector{T}, ispec::Int) where T <: AbstractCompositeBin
    out = getscaled.(bins, ispec)
    i = findall(.!isnothing.(getindex.(out, 1)))
    return (getindex.(out[i], 1), getindex.(out[i], 2))
end


function composite_variance_func(bins::Vector{LinCompositeBin})
    @assert all(nn.(bins) .>= 2)

    Nspec = maximum(maximum(getfield.(bins, :ispec)))
    mm = [Vector{NTuple{2, Int}}() for i in 1:Nspec]
    for j in 1:length(bins)
        i = 1
        for ispec in bins[j].ispec
            push!(mm[ispec], (i, j))
            i += 1
        end
    end

    prog = ProgressUnknown(desc="evaluations:", dt=1.5, showspeed=true, color=:light_black)
    shared = (bins=LogCompositeBin.(bins), prevpars=fill(0., Nspec), Nspec=Nspec, mm=mm, output=fill(0., length(bins)))
    funct = let prog=prog, shared=shared
        params::Vector{Float64} -> begin
            ProgressMeter.next!(prog; showvalues=() -> [(:variance, sum(shared.output .^2))])
            # for j in 1:length(shared.bins)
            #     apply_scale!(shared.bins[j], params)
            #     shared.output[j] = std(shared.bins[j])
            # end
            for ispec in 1:shared.Nspec
                if shared.prevpars[ispec] != params[ispec]
                    for (i, j) in mm[ispec]
                        shared.bins[j].scaled[i] = shared.bins[j].ref[i] + params[ispec]
                    end
                    shared.prevpars[ispec] = params[ispec]
                end
            end
            shared.output .= std.(shared.bins)
            return shared.output
        end
    end
    return prog, shared, funct
end


struct Composite
    R::Union{Nothing, Float64}
    dl::Union{Nothing, Float64}
    refwl::Union{Nothing, Float64}
    bins::Vector{LinCompositeBin}

    function Composite(specs::Vector{SingleSpec};
                       rev=false,      # Reverse redshift ordering
                       R=nothing,      # Resolution
                       dl=nothing,     # Fixed delta-lambda
                       usemodel=false, # Use model rather than data
                       renorm=true,    # Normalize each spectrum before stacking
                       refwl=nothing,  # Use a specific wavelength for normalization
                       fit=false,      # Estimate individual spectra normalization by requiring the output scatter is minimized
                       plot=false)     # Do plot while accumulating spectra
        if usemodel
            @assert  all([all(s.m .> 0) for s in specs])
        else
            @assert  all([all(s.y .> 0) for s in specs])
        end
        logstep(R) = log10(1/R + 1)

        # Sort spectra by redshift
        i = sortperm( getfield.(specs, :z), rev=rev)
        zr = [extrema(getfield.(specs, :z))...]
        specs = specs[i]

        # Prepare domain for the composite spectrum
        xr = [minimum([minimum(s.x) for s in specs]),
              maximum([maximum(s.x) for s in specs])]
        if isnothing(dl)  &&  !isnothing(R)
            domain = 10. .^collect(log10(xr[1]):logstep(R):log10(xr[2]))
        elseif !isnothing(dl)  &&  isnothing(R)
            domain = collect(xr[1]:dl:xr[2])
        else
            error("Only one among R and dl is supposed to be used")
        end

        # Prepare CopositeBin structures
        bins = LinCompositeBin.(domain)

        # Loop through spectra
        @showprogress for ispec in 1:length(specs)
            spec = specs[ispec]

            # Identify range of current spectrum
            @assert issorted(spec.x)
            xr = extrema(spec.x)
            if !isnothing(refwl)
                (xr[1] < refwl < xr[2])  ||  continue
            end
            j = findall(xr[1] .< domain .< xr[2])
            @assert length(j) > 0

            # Resample spectrum
            Y = usemodel  ?  spec.m  :  spec.y
            v = Dierckx.Spline1D(spec.x, Y, k=1, bc="error")(domain[j])
            for i in 1:length(j)
                add_spec!(bins[j[i]], ispec, v[i])
            end

            # Scale spectrum
            if renorm
                if !isnothing(refwl)
                    scale = 1 / Dierckx.Spline1D(spec.x, Y, k=1, bc="error")(refwl)
                elseif fit
                    scale = 1 / mean(Y)
                else
                    if ispec == 1
                        scale = 1 / mean(Y)
                    else
                        tmp = mean.(bins[j])
                        scale = mean(tmp) / mean(Y)
                    end
                end

                apply_scale!.(bins[j], scale, ispec)
            end
        end

        i = findall(nn.(bins) .> 2)
        domain = domain[i]
        bins = bins[i]

        if fit
            save_scaled_as_ref!.(bins)
            scatter = std.(bins)
            @gp :aa hist(scatter)

            config = CMPFit.Config()
            config.ftol   = 1.e-1
            config.xtol   = 1.e-1
            config.gtol   = 1.e-1
            config.covtol = 1.e-1

            prog, shared, funct = composite_variance_func(bins)
            bestfit = CMPFit.cmpfit(funct, fill(0., length(specs)), config=config)
            @info bestfit.elapsed bestfit.orignorm bestfit.bestnorm

            apply_scale!.(bins, Ref(10 .^(bestfit.param)))

            scatter = std.(bins)
            h = hist(scatter)
            @gp :- :aa hist_bins(h) hist_weights(h) "w steps t 'After' lw 3"
        end

        # Global scaling
        save_scaled_as_ref!.(bins)
        apply_scale!.(bins, Ref(fill(1 / mean(mean.(bins)), length(specs))))
        save_scaled_as_ref!.(bins)

        if plot
            data = Dict(:redshift => Float64[],
                        :orig_domain => Float64[],
                        :orig_flux => Float64[],
                        :comp_domain => Float64[],
                        :comp_flux => Float64[],
                        :comp_redshift => Float64[])

            for ispec in 1:length(specs)
                spec = specs[ispec]
                append!(data[:redshift], fill(spec.z, length(spec.x)))
                append!(data[:orig_domain], spec.x)
                append!(data[:orig_flux], usemodel  ?  spec.m  :  spec.y)

                wl, scaled = getscaled(bins, ispec)
                append!(data[:comp_domain], wl)
                append!(data[:comp_flux], scaled)
                append!(data[:comp_redshift], fill(spec.z, length(wl)))
            end

            @gp    :Composite "set grid" xlog=true ylog=true :-
            color = v2argb(:roma, data[:redshift], alpha=0.8, range=[extrema(data[:redshift])...])
            @gp :- :Composite data[:orig_domain] data[:orig_flux] color "w d t 'Data' lc rgb var" :-
            color =  v2argb(:roma, data[:comp_redshift], alpha=0.8, range=[extrema(data[:comp_redshift])...])
            @gp :- :Composite data[:comp_domain] data[:comp_flux] color "w d t 'Resampled and scaled' lc rgb var" :-
            @gp :- :Composite domain       mean.(bins)                    "w l t 'Arith. composite' lc rgb 'black' lw 3" :-
            @gp :- :Composite domain 10 .^ mean.(LogCompositeBin.(bins))  "w l t 'Geom. composite' lc rgb 'black' dt 2 lw 3" :-
            @gp :- :Composite "set autoscale fix"
        end

        return new(R, dl, refwl, bins)
    end
end

domain(  cc::Composite) =  getfield.(cc.bins, :wavelength)
mean(    cc::Composite) =  mean.(    cc.bins)
std(     cc::Composite) =  std.(     cc.bins)
nn(      cc::Composite) =  nn.(      cc.bins)
geommean(cc::Composite) =  mean.(LogCompositeBin.(cc.bins))
geomstd( cc::Composite) =  std.( LogCompositeBin.(cc.bins))



function plot(specs::Vector{SingleSpec}, cc::Composite)
    dom = domain(cc)
    xr = [extrema(dom)...]
    @gp "set grid" "set autoscale noextend" :-
    @gp :- xlog=true xr=xr "set multiplot layout 3,1" "unset colorbox" xlab="" lma=0.06  rma=0.96 :-

    z = getfield.(specs, :z)
    h = hist(z)
    zbins = range([extrema(z)...]..., 200)
    xx = fill(NaN,  length(dom) * length(zbins))
    yy = fill(NaN,  length(dom) * length(zbins))
    zz = fill(NaN, (length(dom),  length(zbins)))
    for i in 1:length(zbins)
        i1 = argmin(abs.(12000 / (zbins[i] + 1) .- dom))
        i2 = argmin(abs.(18000 / (zbins[i] + 1) .- dom))
        zz[i1:i2, i] .= Dierckx.Spline1D(hist_bins(h), hist_weights(h), k=1, bc="error")(zbins[i])
        xx[((i-1) * length(dom) + 1):(i * length(dom))] .= dom
        yy[((i-1) * length(dom) + 1):(i * length(dom))] .= zbins[i]
    end
    @gp :- 1 ylab="Redshift"         "set xtics format ''"       bma=0.78  tma=0.98 yr=extrema(z) ylog=false :-
    @gp :- 1 "set y2tics"            "set ytics nomirror" :-
    @gp :- 1 xx yy zz[:]                                                 "w p            notit                         lc palette pt 4 ps 0.25"  :-
    @gp :- 1 dom nn(cc)                                                  "w l            t 'N. spectra'           lw 3 lc rgb 'black' axes x1y2" :-

    acomp  = mean(cc)
    sacomp = std(cc)
    gcomp  = geommean(cc)
    sgcomp = geomstd(cc)

    @gp :- 2 "set y2label ''"        "set ytics mirror"       "set style fill transparent solid 0.5" :-
    @gp :- 2 ylab="Lum [arb. units]" "set xtics format ''"       bma=0.42  tma=0.77 yr=extrema(acomp) ylog=true :-
    @gp :- 2 dom 10 .^gcomp                       "w l            t 'Geometric composite'  lw 1 lc rgb 'red' dt 3"  :-
    @gp :- 2 dom acomp .- sacomp acomp .+ sacomp  "w filledcurves t 'Arithmetic scatter'        lc rgb 'gray'" :-
    @gp :- 2 dom acomp                                           "w l            t 'Arithmetic composite' lw 1 lc rgb 'blue'" :-

    @gp :- 3 xlab="Wavelength [A]"   "set xtics format '% h'" "set style fill transparent solid 0.5"  :-
    @gp :- 3 ylab="log. Lum [arb. units]"                        bma=0.06  tma=0.41 yr=extrema(gcomp) ylog=false :-
    @gp :- 3 dom log10.(acomp)                    "w l            t 'Arithmetic composite' lw 1 lc rgb 'blue' dt 3" :-
    @gp :- 3 dom gcomp .- sgcomp gcomp .+ sgcomp  "w filledcurves t 'Geometric scatter' lc rgb 'gray'" :-
    @gp :- 3 dom gcomp                            "w l            t 'Geometric composite'  lw 1 lc rgb 'red'"  :-

    cont_x0    = 3000
    cont_norm  = 0.28
    cont_slope = -1.7
    xx = range(xr..., 100)
    @gp :- 3 xx cont_norm .+ cont_slope .* log10.(xx ./ cont_x0) "w l t 'Slope=$(cont_slope)' lw 3 dt 2" :-

    cont_x0    = 6000
    cont_norm  = -0.2
    cont_slope = -1.0
    xx = range(xr..., 100)
    @gp :- 3 xx cont_norm .+ cont_slope .* log10.(xx ./ cont_x0) "w l t 'Slope=$(cont_slope)' lw 3 dt 2"
end


# ====================================================================
serialize_filename = "composite.ser"
if !isfile(serialize_filename)
    input_path  = "input/input_Euclid"
    output_path = "results_Euclid"

    # Read input catalog
    f = FITS("$(input_path)/catalog.fits")
    catalog = DataFrame(f[2])
    close(f)

    specs = Vector{SingleSpec}()
    for i in 1:nrow(catalog)
        try
            filename = "$(output_path)/JSON/$(catalog[i, :object_id]).json.gz"
            res = QSFit.deserialize(filename)
            push!(specs, SingleSpec(coords(res.data.domain),
                                    values(res.data),
                                    res.bestfit(),
                                    catalog[i, :Z]))
        catch err
        end
    end

    # Cut in redshift to avoid useless noise at blue edge
    i = findall(getfield.(specs, :z) .< 4)
    specs = specs[i]

    cc  = Composite(specs, R=2000)
    rcc = Composite(specs, R=2000, rev=true)
    serialize(serialize_filename, (specs, cc, rcc))
else
    specs, cc, rcc = deserialize(serialize_filename)
end


# ====================================================================
aaa()

yy = CSV.read("/home/gcalderone/tmp/Yuming/q1_qsocomp_spec_constant_r500_20251114.csv", DataFrame)
bb = CSV.read("/home/gcalderone/tmp/Yuming/sdss_all_mean_hostcorr.dat", DataFrame);
f = FITS("/home/gcalderone/tmp/Yuming/Salvatore_composite_median_flux_normalization.fits")
ss = DataFrame(f[2])
close(f)
ss = ss[findall(isfinite.(ss.specMean)), :]

refwl = 5600.
@gp    :cmp "set grid" xlabel="Wavelength [A] (rest frame)" ylabel="{/Symbol l} L_{/Symbol l} (arb.units)" xlog=true ylog=true :-
@gp :- :cmp  domain(cc)    domain(cc)     .*       mean(cc)     ./ Dierckx.Spline1D( domain(cc)  ,      mean(cc)     , k=1, bc="error")(refwl) ./ refwl "w l t 'arith'"
@gp :- :cmp domain(rcc)    domain(rcc)    .*      mean(rcc)     ./ Dierckx.Spline1D(domain(rcc)  ,      mean(rcc)    , k=1, bc="error")(refwl) ./ refwl "w l t 'arith rev'"
@gp :- :cmp  domain(cc)    domain(cc)     .* 10 .^geommean( cc) ./ Dierckx.Spline1D( domain(cc)  , 10 .^geommean(cc) , k=1, bc="error")(refwl) ./ refwl "w l t 'geom'"
@gp :- :cmp domain(rcc)    domain(rcc)    .* 10 .^geommean(rcc) ./ Dierckx.Spline1D(domain(rcc)  , 10 .^geommean(rcc), k=1, bc="error")(refwl) ./ refwl "w l t 'geom rev'"
@gp :- :cmp yy.wavelength yy.wavelength   .* yy.mean_flux       ./ Dierckx.Spline1D(yy.wavelength, yy.mean_flux      , k=1, bc="error")(refwl) ./ refwl "w l t 'Yuming (arith)'"
@gp :- :cmp yy.wavelength yy.wavelength   .* yy.geo_flux        ./ Dierckx.Spline1D(yy.wavelength, yy.geo_flux       , k=1, bc="error")(refwl) ./ refwl "w l t 'Yuming (geom)'"
@gp :- :cmp ss.wavelength                    ss.specMean        ./ Dierckx.Spline1D(ss.wavelength, ss.specMean       , k=1, bc="error")(refwl)          "w l t 'Salvatore'"
@gp :- :cmp bb[:, 1]                         bb[:, 2]           ./ Dierckx.Spline1D(     bb[:, 1], bb[:, 2]          , k=1, bc="error")(refwl)          "w l t 'Beta'"

plot(specs, cc)


# Analysis with dedicated recipe
abstract type EuclidComposite <: QSFit.QSORecipes.Type1 end

import QSFit: line_component
function line_component(recipe::CRecipe{<: EuclidComposite}, tid::Union{Float64, Val{TID}}, template::Type{<: QSFit.NarrowLine}) where TID
    @track_recipe
    comp = @invoke line_component(recipe::CRecipe{<: supertype(EuclidComposite)}, tid, template)
    comp.fwhm.high = 3e3
    return comp
end

import QSFit.QSORecipes: add_qso_continuum!
function add_qso_continuum!(recipe::CRecipe{<: EuclidComposite}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    data = QSFit.QSORecipes.getdata(fp, ith)
    λ = coords(GModelFit.domain(data))
    comp = QSFit.sbpl(3000)
    comp.x0.val = median(λ)
    comp.x0.fixed = false
    comp.norm.val = median(values(data)) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    (comp.norm.val < 0)  &&  (comp.norm.val = mad(values(data)))
    @assert comp.norm.val > 0 "Continuum guess normalization is zero or negative"
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha1.val  = -1.5
    comp.alpha1.low  = -5
    comp.alpha1.high =  5
    comp.alpha2.val  = -1.5
    comp.alpha2.low  = -5
    comp.alpha2.high =  5

    QSFit.QSORecipes.getmodel(fp, ith)[:QSOcont] = comp
    push!(QSFit.QSORecipes.getmodel(fp, ith)[:Continuum].list, :QSOcont)
end

# Prepare spectrum
i = findall(nn(cc) .>= 2);
xx = domain(cc)[i]
yy = 10 .^mean(cc, geom=true)[i] .* 1e-17
ee = 10 .^std( cc, geom=true)[i] ./ sqrt.(nn(cc)[i]) .* 1e-17
spec = Spectrum(xx, yy, ee, resolution=1000.)

recipe = CRecipe{EuclidComposite}()
# Consider all range observed by NISP
recipe.wavelength_range = [1215, 2e4]
# Allow nuisance lines only in the region where the iron template is missing
recipe.nuisance_avoid = [[0, 3100], [4000, 1e6]]
# Allow extended exploration of parameter space (i.e. avoid getting stuck when initial values are too far from optimum)
recipe.solver.config.stepfactor = 200

res = analyze(recipe, spec)
viewer(res)
