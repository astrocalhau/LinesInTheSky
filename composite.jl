using Revise, FITSIO, DataFrames, Gnuplot, Dierckx, ProgressMeter, Statistics, StatsBase, Serialization
using QSFit, QSFit.QSORecipes, GModelFit, GModelFitViewer

Gnuplot.options.term = "qt size 1800,1300 enhanced font 'Latin Modern Roman, 13' lw 1.5"

struct SingleSpec
    x::Vector{Float64}
    y::Vector{Float64}
    m::Vector{Float64}
    z::Float64
end


collapse(m; kws...) = collapse(m, 1:size(m)[2]; kws...)
function collapse(m, j; robust=false)
    avg = fill(NaN, length(j))
    sig = fill(NaN, length(j))
    nn  = fill(0 , length(j))
    c = 1
    for j in j
        i = findall(.!isnan.(m[:, j]))
        if length(i) > 0
            nn[c] = length(i)
            if robust
                avg[c] = median(m[i, j])
                sig[c] = mad(   m[i, j])
            else
                avg[c] = mean(  m[i, j])
                sig[c] = std(   m[i, j])
            end
        end
        c += 1
    end
    return avg, sig, nn
end

struct Composite
    domain::Vector{Float64}
    composite::Vector{Float64}
    scatter::Vector{Float64}
    geom_composite::Vector{Float64}
    geom_scatter::Vector{Float64}
    nn::Vector{Float64}

    function Composite(specs::Vector{SingleSpec};
                       rev=false,      # Reverse redshift ordering
                       R=nothing,      # Resolution
                       dl=nothing,     # Fixed delta-lambda
                       usemodel=false, # Use model rather than data
                       renorm=true,    # Normalize each spectrum before stacking
                       refwl=nothing,  # Use a specific wavelength for normalization
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

        # Prepare domain for the composite spectrum and matrix to store all rebinned and scaled spectra
        xr = [minimum([minimum(s.x) for s in specs]),
              maximum([maximum(s.x) for s in specs])]
        if isnothing(dl)  &&  !isnothing(R)
            domain = 10. .^collect(log10(xr[1]):logstep(R):log10(xr[2]))
        elseif !isnothing(dl)  &&  isnothing(R)
            domain = collect(xr[1]:dl:xr[2])
        else
            error("Only one among R and dl is supposed to be used")
        end
        matrix = fill(NaN, (length(specs), length(domain)))
        plot_data = Dict(:redshift => Float64[],
                         :orig_domain => Float64[],
                         :orig_flux => Float64[],
                         :comp_domain => Float64[],
                         :comp_flux => Float64[],
                         :comp_redshift => Float64[])

        # Loop through spectra
        @showprogress for i in 1:length(specs)
            spec = specs[i]

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
            resampledY = Dierckx.Spline1D(spec.x, Y, k=1, bc="error")(domain[j])

            # Scale spectrum
            if renorm
                if !isnothing(refwl)
                    scaledY = resampledY ./ Dierckx.Spline1D(spec.x, Y, k=1, bc="error")(refwl)
                else
                    if i == 1
                        scale = 1.0
                    else
                        tmp, _, _ = collapse(matrix, j)
                        scale = mean(tmp[findall(.!isnan.(tmp))])
                    end

                    scaledY = scale .* resampledY ./ mean(Y)
                end
            else
                scaledY = resampledY
            end

            # Store resampled and scaled spectrum
            matrix[i, j] .= scaledY

            if plot
                append!(plot_data[:redshift], fill(spec.z, length(spec.x)))
                append!(plot_data[:orig_domain], spec.x)
                append!(plot_data[:orig_flux], Y)
                append!(plot_data[:comp_domain], domain[j])
                append!(plot_data[:comp_flux], scaledY)
                append!(plot_data[:comp_redshift], fill(spec.z, length(j)))
                # @gp :aa "set autoscale fix" xlog=true ylog=true domain[j] scaledY "w p" spec.x Y "w p"
                # readline()
            end
        end

        composite          , scatter, nn = collapse(matrix)
        geom_composite, geom_scatter, _  = collapse(log10.(matrix))

        scale = mean(composite)
        plot_data[:comp_flux] ./= scale
        composite             ./= scale
        scatter               ./= scale
        geom_composite        .-= log10(scale)  # no need to apply offset on geom_scatter here

        if plot
            @gp "set grid" xlog=true ylog=true :-
            color = v2argb(:roma, plot_data[:redshift], alpha=0.8, range=[extrema(plot_data[:redshift])...])
            @gp :- plot_data[:orig_domain] plot_data[:orig_flux] color "w d t 'Data' lc rgb var" :-
            color = v2argb(:roma, plot_data[:comp_redshift], alpha=0.8, range=[extrema(plot_data[:comp_redshift])...])
            @gp :- plot_data[:comp_domain] plot_data[:comp_flux] color "w d t 'Resampled and scaled' lc rgb var" :-
            @gp :- domain composite            "w l t 'Arith. composite' lc rgb 'black' lw 3" :-
            @gp :- domain 10 .^ geom_composite "w l t 'Geom. composite' lc rgb 'black' dt 2 lw 3" :-
            @gp :- "set autoscale fix"
        end
        return new(domain, composite, scatter, geom_composite, geom_scatter, nn)
    end
end



function plot(specs::Vector{SingleSpec}, cc::Composite)
    xr = [extrema(cc.domain)...]
    @gp "set grid" "set autoscale noextend" :-
    @gp :- xlog=true xr=xr "set multiplot layout 3,1" "unset colorbox" xlab="" lma=0.06  rma=0.96 :-

    z = getfield.(specs, :z)
    h = hist(z)
    zbins = range([extrema(z)...]..., 200)
    xx = fill(NaN,  length(cc.domain) * length(zbins))
    yy = fill(NaN,  length(cc.domain) * length(zbins))
    zz = fill(NaN, (length(cc.domain),  length(zbins)))
    for i in 1:length(zbins)
        i1 = argmin(abs.(12000 / (zbins[i] + 1) .- cc.domain))
        i2 = argmin(abs.(18000 / (zbins[i] + 1) .- cc.domain))
        zz[i1:i2, i] .= Dierckx.Spline1D(hist_bins(h), hist_weights(h), k=1, bc="error")(zbins[i])
        xx[((i-1) * length(cc.domain) + 1):(i * length(cc.domain))] .= cc.domain
        yy[((i-1) * length(cc.domain) + 1):(i * length(cc.domain))] .= zbins[i]
    end
    @gp :- 1 ylab="Redshift"         "set xtics format ''"       bma=0.78  tma=0.98 yr=extrema(z) ylog=false :-
    @gp :- 1 "set y2tics"            "set ytics nomirror" :-
    @gp :- 1 xx yy zz[:]                                                      "w p            notit                         lc palette pt 4 ps 0.25"  :-
    @gp :- 1 cc.domain cc.nn                                                  "w l            t 'N. spectra'           lw 3 lc rgb 'black' axes x1y2" :-

    @gp :- 2 "set y2label ''"        "set ytics mirror"       "set style fill transparent solid 0.5" :-
    @gp :- 2 ylab="Lum [arb. units]" "set xtics format ''"       bma=0.42  tma=0.77 yr=extrema(cc.composite) ylog=true :-
    @gp :- 2 cc.domain 10 .^cc.geom_composite                                 "w l            t 'Geometric composite'  lw 1 lc rgb 'red' dt 3"  :-
    @gp :- 2 cc.domain cc.composite .- cc.scatter cc.composite .+ cc.scatter  "w filledcurves t 'Arithmetic scatter'        lc rgb 'gray'" :-
    @gp :- 2 cc.domain cc.composite                                           "w l            t 'Arithmetic composite' lw 1 lc rgb 'blue'" :-

    @gp :- 3 xlab="Wavelength [A]"   "set xtics format '% h'" "set style fill transparent solid 0.5"  :-
    @gp :- 3 ylab="log. Lum [arb. units]"                        bma=0.06  tma=0.41 yr=extrema(cc.geom_composite) ylog=false :-
    @gp :- 3 cc.domain log10.(cc.composite)                                   "w l            t 'Arithmetic composite' lw 1 lc rgb 'blue' dt 3" :-
    @gp :- 3 cc.domain cc.geom_composite .- cc.geom_scatter cc.geom_composite .+ cc.geom_scatter  "w filledcurves t 'Geometric scatter' lc rgb 'gray'" :-
    @gp :- 3 cc.domain cc.geom_composite                                      "w l            t 'Geometric composite'  lw 1 lc rgb 'red'"  :-

    cont_x0    = 3000
    cont_norm  = 0.28
    cont_slope = -1.65
    xx = range(xr..., 100)
    @gp :- 3 xx cont_norm .+ cont_slope .* log10.(xx ./ cont_x0) "w l t 'Slope=$(cont_slope)' lw 3 dt 2" :-

    cont_x0    = 6000
    cont_norm  = -0.18
    cont_slope = -1.05
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



yy = CSV.read("/home/gcalderone/tmp/Yuming/q1_qsocomp_spec_constant_r500_20251114.csv", DataFrame)
bb = CSV.read("/home/gcalderone/tmp/Yuming/sdss_all_mean_hostcorr.dat", DataFrame);
refwl = 5600.
@gp    :cmp "set grid" xlog=true ylog=true :-
@gp :- :cmp  cc.domain     cc.domain    .*  cc.composite           ./ Dierckx.Spline1D( cc.domain    ,  cc.composite          , k=1, bc="error")(refwl) ./ refwl "w l t 'arith'"
@gp :- :cmp rcc.domain    rcc.domain    .* rcc.composite           ./ Dierckx.Spline1D(rcc.domain    , rcc.composite          , k=1, bc="error")(refwl) ./ refwl "w l t 'arith rev'"
@gp :- :cmp  cc.domain     cc.domain    .* 10 .^ cc.geom_composite ./ Dierckx.Spline1D( cc.domain    , 10 .^ cc.geom_composite, k=1, bc="error")(refwl) ./ refwl "w l t 'geom'"
@gp :- :cmp rcc.domain    rcc.domain    .* 10 .^rcc.geom_composite ./ Dierckx.Spline1D(rcc.domain    , 10 .^rcc.geom_composite, k=1, bc="error")(refwl) ./ refwl "w l t 'geom rev'"
@gp :- :cmp yy.wavelength yy.wavelength .* yy.mean_flux            ./ Dierckx.Spline1D(yy.wavelength , yy.mean_flux           , k=1, bc="error")(refwl) ./ refwl "w l t 'Yuming'"
@gp :- :cmp yy.wavelength yy.wavelength .* yy.geo_flux             ./ Dierckx.Spline1D(yy.wavelength , yy.geo_flux            , k=1, bc="error")(refwl) ./ refwl "w l t 'Yuming (geom)'"
@gp :- :cmp bb[:, 1]           bb[:, 1] .* bb[:, 2]                ./ Dierckx.Spline1D(     bb[:, 1] , bb[:, 2]               , k=1, bc="error")(refwl) ./ refwl "w l t 'Beta'"


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
    λ = coords(domain(data))
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
i = findall(cc.nn .>= 2);
xx = cc.domain[i]
yy = 10 .^cc.geom_composite[i] .* 1e-17
ee = 10 .^cc.geom_scatter[i] ./ sqrt.(cc.nn[i]) .* 1e-17
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
