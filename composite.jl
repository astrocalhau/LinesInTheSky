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
    avg = fill(0., length(j))
    sig = fill(0., length(j))
    nn  = fill(0 , length(j))
    c = 1
    for j in j
        i = findall(.!isnan.(m[:, j]))
        nn[c] = length(i)
        if robust
            avg[c] = median(m[i, j])
            sig[c] = mad(   m[i, j])
        else
            avg[c] = mean(m[i, j])
            sig[c] = std( m[i, j])
        end
        c += 1
    end
    return avg, sig, nn
end

struct Composite
    domain::Vector{Float64}
    geom::Bool
    composite::Vector{Float64}
    scatter::Vector{Float64}
    nn::Vector{Float64}

    function Composite(specs::Vector{SingleSpec}; R=nothing, dl=nothing, geom=false, rev=false, usemodel=false)
        if !geom
            if usemodel
                @assert  all([all(s.m .> 0) for s in specs])
            else
                @assert  all([all(s.y .> 0) for s in specs])
            end
        end
        logstep(R) = log10(1/R + 1)

        # Sort spectra by redshift
        i = sortperm(getfield.(specs, :z), rev=rev)
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

        # Loop through spectra
        @showprogress for i in 1:length(specs)
            spec = specs[i]

            # Identify range of current spectrum
            xr = extrema(spec.x)
            j = findall(xr[1] .< domain .< xr[2])
            @assert length(j) > 0

            Y = usemodel  ?  spec.m  :  spec.y
            if !geom
                # Scale spectrum
                if i == 1
                    scale = 1.0
                else
                    tmp, _, _ = collapse(matrix, j)
                    scale = mean(tmp[findall(.!isnan.(tmp))])
                end
                scale /= mean(Y)

                # Resample and store
                matrix[i, j] .= scale .* Dierckx.Spline1D(spec.x, Y, k=1, bc="error")(domain[j])
            else
                # Scale spectrum
                if i == 1
                    off = 0.0
                else
                    tmp, _, _ = collapse(matrix, j)
                    off = mean(tmp[findall(.!isnan.(tmp))])
                end
                off -= mean(log10.(Y))

                # Resample and store
                matrix[i, j] .= off .+ Dierckx.Spline1D(spec.x, log10.(Y), k=1, bc="error")(domain[j])
            end
        end

        if !geom
            matrix ./= mean(matrix[.!isnan.(matrix)])
        else
            matrix .-= mean(matrix[.!isnan.(matrix)])
        end

        composite, scatter, nn = collapse(matrix, robust=false)

        if !geom
            return new(domain, geom, composite, scatter, nn)
        else
            return new(domain, geom, 10 .^composite, 10 .^scatter, nn)
        end
    end
end



function plot(specs::Vector{SingleSpec}, ac::Composite, gc::Composite)
    @assert !ac.geom
    @assert gc.geom

    xr = [extrema(ac.domain)...]
    @gp "set grid" "set autoscale noextend" :-
    @gp :- xlog=true xr=xr "set multiplot layout 3,1" "unset colorbox" xlab="" lma=0.06  rma=0.96 :-

    z = getfield.(specs, :z)
    h = hist(z)
    zbins = range([extrema(z)...]..., 200)
    xx = fill(NaN,  length(ac.domain) * length(zbins))
    yy = fill(NaN,  length(ac.domain) * length(zbins))
    zz = fill(NaN, (length(ac.domain),  length(zbins)))
    for i in 1:length(zbins)
        i1 = argmin(abs.(12000 / (zbins[i] + 1) .- ac.domain))
        i2 = argmin(abs.(18000 / (zbins[i] + 1) .- ac.domain))
        zz[i1:i2, i] .= Dierckx.Spline1D(hist_bins(h), hist_weights(h), k=1, bc="error")(zbins[i])
        xx[((i-1) * length(ac.domain) + 1):(i * length(ac.domain))] .= ac.domain
        yy[((i-1) * length(ac.domain) + 1):(i * length(ac.domain))] .= zbins[i]
    end
    @gp :- 1 ylab="Redshift"         "set xtics format ''"       bma=0.78  tma=0.98 yr=extrema(z) ylog=false :-
    @gp :- 1 "set y2tics"            "set ytics nomirror" :-
    @gp :- 1 xx yy zz[:]                                                      "w p            notit                         lc palette pt 4 ps 0.25"  :-
    @gp :- 1 ac.domain ac.nn                                                  "w l            t 'N. spectra'           lw 3 lc rgb 'black' axes x1y2" :-

    @gp :- 2 "set y2label ''"        "set ytics mirror"       "set style fill transparent solid 0.5" :-
    @gp :- 2 ylab="Lum [arb. units]" "set xtics format ''"       bma=0.42  tma=0.77 yr=extrema(ac.composite) ylog=true :-
    @gp :- 2 gc.domain gc.composite                                           "w l            t 'Geometric composite'  lw 1 lc rgb 'red' dt 3"  :-
    @gp :- 2 ac.domain ac.composite .- ac.scatter ac.composite .+ ac.scatter  "w filledcurves t 'Arithmetic scatter'        lc rgb 'gray'" :-
    @gp :- 2 ac.domain ac.composite                                           "w l            t 'Arithmetic composite' lw 1 lc rgb 'blue'" :-

    @gp :- 3 xlab="Wavelength [A]"   "set xtics format '% h'" "set style fill transparent solid 0.5"  :-
    @gp :- 3 ylab="log. Lum [arb. units]"                        bma=0.06  tma=0.41 yr=extrema(gc.composite) ylog=true :-
    @gp :- 3 ac.domain ac.composite                                           "w l            t 'Arithmetic composite' lw 1 lc rgb 'blue' dt 3" :-
    @gp :- 3 gc.domain gc.composite ./ gc.scatter gc.composite .* gc.scatter  "w filledcurves t 'Geometric scatter'         lc rgb 'gray'" :-
    @gp :- 3 gc.domain gc.composite                                           "w l            t 'Geometric composite'  lw 1 lc rgb 'red'"  :-

    cont_x0    = 3000
    cont_norm  = 10^0.4
    cont_slope = -1.6
    xx = range(xr..., 100)
    @gp :- 3 xx cont_norm .* (xx ./ cont_x0).^cont_slope "w l t 'Slope=$(cont_slope)' lw 3 dt 2" :-

    cont_x0    = 6000
    cont_norm  = 0.88
    cont_slope = -1.0
    xx = range(xr..., 100)
    @gp :- 3 xx cont_norm .* (xx ./ cont_x0).^cont_slope "w l t 'Slope=$(cont_slope)' lw 3 dt 2"
end


# Euclid
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

    ac  = Composite(specs, R=2000)
    gc  = Composite(specs, R=2000, geom=true)
    rac = Composite(specs, R=2000           , rev=true)
    rgc = Composite(specs, R=2000, geom=true, rev=true)

    serialize(serialize_filename, (specs, ac, gc, rac, rgc))
else
    specs, ac, gc, rac, rgc = deserialize(serialize_filename)
end
    
@gp    :cmp "set grid" ac.domain ac.composite "w l t 'arith'" rac.domain rac.composite "w l t 'arith rev'" xlog=true ylog=true :-
@gp :- :cmp            gc.domain gc.composite "w l t 'geom'"  rgc.domain rgc.composite "w l t 'geom rev'"
plot(specs, ac, gc)



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
i = findall(gc.nn .>= 2);
spec = Spectrum(gc.domain[i], gc.composite[i] .* 1e-17, gc.scatter[i] ./ sqrt.(gc.nn[i]) .* 1e-17, resolution=1000.)

recipe = CRecipe{EuclidComposite}()
# Consider all range observed by NISP
recipe.wavelength_range = [1215, 2e4]
# Allow nuisance lines only in the region where the iron template is missing
recipe.nuisance_avoid = [[0, 3100], [4000, 1e6]]
# Allow extended exploration of parameter space (i.e. avoid getting stuck when initial values are too far from optimum)
recipe.solver.config.stepfactor = 200

res = analyze(recipe, spec)
viewer(res)

