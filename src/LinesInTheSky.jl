module LinesInTheSky

using GModelFit, QSFit, QSFit.ATL, QSFit.QSORecipes, Statistics, StatsBase

export WP9Type1IR



import QSFit: init_recipe!, add_line!
import QSFit.ATL: Permitted, Forbidden
import QSFit.QSORecipes: lines_dict
using  QSFit.QSORecipes:  getmodel, getdomain, getdata


function __init__()
    ATL.register(Permitted,  "HeIλ10832"  , [(10832.0574, 19.819635, 20.964240),    # He I      |E1|                        1s.2s-1s.2p                        |       3S-3Po        |   1-0   |   19.819635 -   20.964240|
                                             (10833.2167, 19.819635, 20.964117),    # He I      |E1|                        1s.2s-1s.2p                        |       3S-3Po        |   1-1   |   19.819635 -   20.964117|
                                             (10833.3064, 19.819635, 20.964108)])   # He I      |E1|                        1s.2s-1s.2p                        |       3S-3Po        |   1-2   |   19.819635 -   20.964108|

    ATL.register(Forbidden,  "[SIII]λ9069",  ( 9071.1   ,  0.037031, 1.403836))     # [S III]   |M1|                       3s2.3p2-3s2.3p2                     |          3P-1D      |   1-2   |     0.037031 -   1.403836|
    ATL.register(Forbidden,  "[SIII]λ9532",  ( 9533.2   ,  0.103287, 1.403836))     # [S III]   |M1|                       3s2.3p2-3s2.3p2                     |          3P-1D      |   2-2   |     0.103287 -   1.403836|

    ATL.register(Permitted,  "Pa9"        ,  ( 9231.547 , 12.087510, 13.430559))    # H I       |E1|                           3*-9*                           |         3-9         |   *-*   |   12.087510 -   13.430559|
    ATL.register(Permitted,  "Pa8"        ,  ( 9548.590 , 12.087510, 13.385966))    # H I       |E1|                           3*-8*                           |         3-8         |   *-*   |   12.087510 -   13.385966|
    ATL.register(Permitted,  "Paδ"        ,  (10052.128 , 12.087510, 13.320923))    # H I       |E1|                           3*-7*                           |         3-7         |   *-*   |   12.087510 -   13.320923|
    ATL.register(Permitted,  "Paγ"        ,  (10941.091 , 12.087510, 13.220709))    # H I       |E1|                           3*-6*                           |         3-6         |   *-*   |   12.087510 -   13.220709|
    ATL.register(Permitted,  "Paβ"        ,  (12821.59  , 12.087510, 13.054507))    # H I       |E1|                           3*-5*                           |         3-5         |   *-*   |   12.087510 -   13.054507|
    ATL.register(Permitted,  "Paα"        ,  (18756.13  , 12.087510, 12.748544))    # H I       |E1|                           3*-4*                           |         3-4         |   *-*   |   12.087510 -   12.748544|
    ATL.register(Permitted,  "OIλ8448"    , [(8448.568  ,  9.521367, 10.988885),    # O I		|E1|	2s2.2p3.(4So).3s-2s2.2p3.(4So).3p	|	3So-3P	|	1-0	|	9.521367 - 10.988885	|
	                                         (8448.680  ,  9.521367, 10.988866),    # O I		|E1|	2s2.2p3.(4So).3s-2s2.2p3.(4So).3p	|	3So-3P	|	1-2	|	9.521367 - 10.988866	|
	                                         (8449.079  ,  9.521367, 10.988796)])   # O I		|E1|	2s2.2p3.(4So).3s-2s2.2p3.(4So).3p	|	3So-3P	|	1-1	|	9.521367 - 10.988796	|

    ATL.register(Permitted,  "OIλ11289"   , [(11289.407 , 10.988796, 12.087032),    # O I		|E1|	2s2.2p3.(4So).3p-2s2.2p3.(4So).3d	|	3P-3Do	|	1-2	|	10.988796 - 12.087032	|
	                                         (11289.496 , 10.988796, 12.087023),    # O I		|E1|	2s2.2p3.(4So).3p-2s2.2p3.(4So).3d	|	3P-3Do	|	1-1	|	10.988796 - 12.087023	|
	                                         (11290.005 , 10.988866, 12.087043),    # O I		|E1|	2s2.2p3.(4So).3p-2s2.2p3.(4So).3d	|	3P-3Do	|	2-3	|	10.988866 - 12.087043	|
	                                         (11290.120 , 10.988866, 12.087032),    # O I		|E1|	2s2.2p3.(4So).3p-2s2.2p3.(4So).3d	|	3P-3Do	|	2-2	|	10.988866 - 12.087032	|
	                                         (11290.209 , 10.988866, 12.087023),    # O I		|E1|	2s2.2p3.(4So).3p-2s2.2p3.(4So).3d	|	3P-3Do	|	2-1	|	10.988866 - 12.087023	|
	                                         (11290.409 , 10.988885, 12.087023)])   # O I		|E1|	2s2.2p3.(4So).3p-2s2.2p3.(4So).3d	|	3P-3Do	|	0-1	|	10.988885 - 12.087023	|
    ATL.register(Permitted, "OIλ13168"    , [(13167.490 , 10.988796, 11.930390),    # O I        E1              2s2.2p3.(4So).3p-2s2.2p3.(4So).4s                      3P-3So           1-1       10.988796 -   11.930390
                                             (13168.459 , 10.988866, 11.930390),    # O I        E1              2s2.2p3.(4So).3p-2s2.2p3.(4So).4s                      3P-3So           2-1       10.988866 -   11.930390
                                             (13168.731 , 10.988885, 11.930390)])   # O I        E1              2s2.2p3.(4So).3p-2s2.2p3.(4So).4s                      3P-3So           0-1       10.988885 -   11.930390

    ATL.register(Permitted, "CaIIλ8500"   ,  ( 8500.36  ,  13650.19,25414.40))      # Ca II   |   E1          |         3p6.(1S).3d-3p6.(1S).4p      |        2D-2Po  |    3/2-3/2    |  13650.19 -    25414.40   |
    ATL.register(Permitted, "CaIIλ8544"   ,  ( 8544.44  ,  13710.88,25414.40))      # Ca II   |   E1          |         3p6.(1S).3d-3p6.(1S).4p      |        2D-2Po  |    5/2-3/2    |  13710.88 -    25414.40   |
    ATL.register(Permitted, "CaIIλ8664"   ,  ( 8664.52  ,  13650.19,25191.51))      # Ca II   |   E1          |         3p6.(1S).3d-3p6.(1S).4p      |        2D-2Po  |    3/2-1/2    |  13650.19 -    25191.51   |

    ATL.register(Permitted, "HIλ9231"     ,  ( 9231.547 ,  97492.31,108324.73))     # H I     |   E1          |                  3*-9*               |        3-9     |       *-*     |   97492.31 -   108324.73  | + Broad component
    ATL.register(Permitted, "HIλ10052"    ,  (10052.128 ,  97492.31,107440.45))     # H I   |   E1          |                  3*-7*               |        3-7     |       *-*     |   97492.31 -   107440.45  | + Broad component
    ATL.register(Permitted, "HIλ10941"    ,  (10941.091 ,  97492.31,106632.17))     # H I   |   E1          |                  3*-6*               |        3-6     |       *-*     |   97492.31 -   106632.17  | + Broad component
    ATL.register(Permitted, "HIλ12821"    ,  (12821.59  ,  97492.31,105291.66))     #  H I   |   E1          |                  3*-5*               |        3-5     |       *-*     |   97492.31 -   105291.66  | + Broad component
    ATL.register(Permitted, "HIλ18756"    ,  (18756.13  ,  97492.31,102823.90))     #  H I   |   E1          |                  3*-4*               |        3-4     |       *-*     |   97492.31 -   102823.90  | + Broad component
    ATL.register(Permitted, "HIλ19450"    ,  (19450.871 ,  102823.90,107965.06))    # H I  |   E1          |                  4*-8*               |        4-8     |       *-*     |  102823.90 -   107965.06  | + Broad component

    ATL.register(Forbidden, "SVIIIλ9913"  ,  ( 9913.8   ,  0.00,   10087.00))       #   [S VIII] |    M1      |                 2s2.2p5-2s2.2p5     |                         2Po-2Po     |    3/2-1/2      |    0.00 -    10087.00   |

    ATL.register(Forbidden, "SIXλ12523"   ,  (12523.0   ,  0.00,   7985.00))        #    [S IX]   |   M1        |               2s2.2p4-2s2.2p4       |                        3P-3P       |     2-1      |      0.00 -     7985.00    |

    ATL.register(Forbidden, "CIλ9826"     ,  ( 9826.82  , 16.40,  10192.63))        #       [C I]   |    M1           |            2s2.2p2-2s2.2p2       |                        3P-1D      |      1-2      |     16.40 -    10192.63   |

    ATL.register(Forbidden, "CIλ9852"     ,  ( 9852.96  , 43.40,  10192.63))        #       [C I]    |   M1             |          2s2.2p2-2s2.2p2        |                      3P-1D      |      2-2      |     43.40 -    10192.63   |

    ATL.register(Permitted, "HeIIλ10126"  ,  (10126.39  , 51.016783,52.241153))     #       He II   |   E1               |             4*-5*        |                  4-5      |      *-*    |   51.016783 -   52.241153   | + Broad
    ATL.register(Permitted, "HeIIλ11629"  ,  (11629.60  , 52.241153,53.307263))     #        He II      E1                            5*-7*                                     5-7             *-*       52.241153 -   53.307263

    ATL.register(Forbidden, "SIIλ10289"   ,  (10289.55  ,  1.841531,3.046484))      #       [S II]   |   M1             |          3s2.3p3-3s2.3p3       |                2Do-2Po    |     3/2-3/2   |   1.841531 -    3.046484   |
    ATL.register(Forbidden, "SIIλ10323"   ,  (10323.32  ,  1.845472,3.046484))      #       [S II]   |   M1            |           3s2.3p3-3s2.3p3            |                  2Do-2Po     |    5/2-3/2   |   1.845472 -    3.046484   |
    ATL.register(Forbidden, "SIIλ10339"   ,  (10339.24  ,  1.841531,3.040693))      #       [S II]   |   M1            |           3s2.3p3-3s2.3p3            |                  2Do-2Po     |    3/2-1/2   |   1.841531 -    3.040693   |
    ATL.register(Forbidden, "SIIλ10373"   ,  (10373.34  ,  1.845472,3.040693))      #       [S II]   |   E2            |           3s2.3p3-3s2.3p3       |                    2Do-2Po    |     5/2-1/2   |   1.845472 -    3.040693   |

    ATL.register(Forbidden, "NIλ10400"    ,  (10400.587 ,  2.383531,3.575620))      #    [N I]       M1                       2s2.2p3-2s2.2p3                              2Do-2Po         5/2-3/2      2.383531 -    3.575620
    ATL.register(Forbidden, "NIλ10410"    , [(10410.021 ,  2.384611,3.575620),      #      [N I]   |    M1             |          2s2.2p3-2s2.2p3             |                 2Do-2Po    |     3/2-3/2   |   2.384611 -    3.575620   |
	                                         (10410.439 ,  2.384611,3.575572)])     #      [N I]    |   M1            |           2s2.2p3-2s2.2p3           |                   2Do-2Po    |     3/2-1/2  |   2.384611 -    3.575572   |

    ATL.register(Permitted, "HeIλ10832"   , [(10832.0574, 19.819635,20.964240),     #      He I   |    E1           |            1s.2s-1s.2p             |                   3S-3Po     |     1-0    |   19.819635 -   20.964240   |
	                                         (10833.2167, 19.819635,20.964117),     #      He I   |    E1             |            1s.2s-1s.2p           |                      3S-3Po      |     1-1    |   19.819635 -   20.964117   |
	                                         (10833.3064, 19.819635,20.964108)])    #      He I    |   E1            |            1s.2s-1s.2p          |                     3S-3Po       |   1-2   |  19.819635 -   20.964108   | + Broad component

    ATL.register(Forbidden, "PIIλ11471"   ,  (11471.30  , 0.020445,1.101266))       #       [P II]   |   M1             |        3s2.3p2-3s2.3p2        |                     3P-1D       |   1-2   |  0.020445 -    1.101266   |
    ATL.register(Forbidden, "PIIλ11886"   ,  (11886.10  , 0.058163,1.101266))       #       [P II]   |   M1              |         3s2.3p2-3s2.3p2       |                        3P-1D       |    2-2   |  0.058163 -    1.101266   |

    ATL.register(Forbidden, "SiXλ14304"   ,  (14304.9   , 0.000000,0.866724))       #        [Si X]   |   M1             |          2s2.2p-2s2.2p        |                     2Po-2Po      |   1/2-3/2   | 0.000000 -    0.866724   |

    ATL.register(Forbidden, "SiVIλ19650"  ,  (19650.0   , 0.000000,0.631080))       #         [Si VI]  |   M1           |            2s2.2p5-2s2.2p5      |                  2Po-2Po      |   3/2-1/2   |   0.000000 -    0.631080   |
    ATL.register(Permitted, "CIVλ1550", (1550.777, 0.00000, 7.994975)) #            |           C IV   |   E1           |                2S-2Po         1/2-1/2      0.000000 -    7.994975      |

#=
    #Iron lines and multiplets. Need to figure out how to implement them into the iron routines of Giorgio, so they show up on the iron plots.
    ATL.register(Forbidden, "FeIIλ9229"   ,  ( 9229.149 , 0.351861,1.695260))       #      [Fe II]     E2                           3d7-3d7                                  a4F-a4P         5/2-3/2      0.351861 -    1.695260
    ATL.register(Forbidden, "FeIIλ12570"  , [(12570.238 , 0.000000,0.986332),       #      [Fe II]     M1                   3d6.(5D).4s-3d6.(5D).4s                          a6D-a4D         9/2-7/2      0.000000 -    0.986332
                                             (12572.088 , 2.828125,3.814311)])      #      [Fe II]     M1                  3d6.(3F4).4s-3d6.(3G).4s                          b4F-b2G         7/2-7/2      2.828125 -    3.814311
    ATL.register(Forbidden, "FeIIλ12791"  ,  (12791.255 , 0.106950,1.076240))       #      [Fe II]     M1                   3d6.(5D).4s-3d6.(5D).4s                          a6D-a4D         3/2-3/2      0.106950 -    1.076240
    ATL.register(Forbidden, "FeIIλ13209"  ,  (13209.151 , 0.047708,0.986332))       #      [Fe II]     M1                   3d6.(5D).4s-3d6.(5D).4s                          a6D-a4D         7/2-7/2      0.047708 -    0.986332
    ATL.register(Forbidden, "FeIIλ15338"  ,  (15338.903 , 0.232169,1.040468))       #      [Fe II]     E2                           3d7-3d6.(5D).4s                          a4F-a4D         9/2-5/2      0.232169 -    1.040468
    ATL.register(Forbidden, "FeIIλ16439"  ,  (16439.981 , 0.232169,0.986332))       #      [Fe II]     M1                           3d7-3d6.(5D).4s                          a4F-a4D         9/2-7/2      0.232169 -    0.986332
    ATL.register(Forbidden, "FeIIλ16773"  ,  (16773.342 , 0.301294,1.040468))       #      [Fe II]     M1                           3d7-3d6.(5D).4s                          a4F-a4D         7/2-5/2      0.301294 -    1.040468

    #Iron multiplets that might join with other lines
    ATL.register(Permitted, "FeIIλ8929"   ,  ( 8929.086 , 9.849242,11.237786))      #      Fe II      E1                   3d6.(5D).5s-3d6.(5D).5p                          e4D-4Do         7/2-5/2      9.849242 -   11.237786
    ATL.register(Permitted, "FeIIλ9125"   ,  ( 9125.419 , 9.849242,11.207911))      #       Fe II      E1                   3d6.(5D).5s-3d6.(5D).5p                          e4D-4Do         7/2-7/2      9.849242 -   11.207911
    ATL.register(Permitted, "FeIIλ9134"   , [( 9134.54  , 7.772557,9.129870),       #        Fe II      E1                  3d6.(3P4).4p-3d5.4s2                             y4Do-4F          5/2-3/2      7.772557 -    9.129870
                                             ( 9134.85  ,10.544687,11.901953),      #        Fe II      E1                   3d6.(5D).4d-3d5.(4P).4s.4p.(3Po)                 f4D-4Do         3/2-5/2     10.544687 -   11.901953
                                             ( 9134.872 , 9.849242,11.206505)])     #       Fe II      E1                   3d6.(5D).5s-3d6.(5D).5p                          e4D-4Fo         7/2-9/2      9.849242 -   11.206505
    ATL.register(Permitted, "FeIIλ9178"   ,  ( 9178.414 , 9.904543,11.255367))      #       Fe II      E1                   3d6.(5D).5s-3d6.(5D).5p                          e4D-4Fo         5/2-7/2      9.904543 -   11.255367
    ATL.register(Permitted, "FeIIλ9180"   , [( 9180.324 ,11.434514,12.785057),      #       Fe II      E1                   3d6.(5D).5p-3d6.(5D).6s                          4Po-4D          3/2-3/2     11.434514 -   12.785057
                                             ( 9180.568 , 9.940806,11.291313)])     #       Fe II      E1                   3d6.(5D).5s-3d6.(5D).5p                          e4D-4Fo         3/2-5/2      9.940806 -   11.291313
    ATL.register(Permitted, "FeIIλ9199"   ,  (9199.394  , 9.940806,11.288549))      #       Fe II      E1                   3d6.(5D).5s-3d6.(5D).5p                          e4D-4Do         3/2-3/2      9.940806 -   11.288549
    =#
end

abstract type WP9Type1IR <: Type1 end

#to change the allowed max of fwhm
import QSFit.line_component
function line_component(recipe::CRecipe{T}, tid::Union{Val{TID}, Float64}, template::Type{<: ForbiddenLine}) where {TID,T <: WP9Type1IR} 
    @track_recipe
    comp=@invoke line_component(recipe::CRecipe{<: Type1},tid,template)
    comp.fwhm.high=1000.
    comp.fwhm.low=100.
    return comp
end
 

# To use different QSO Continuum models. Can be deleted/commented out if using the standard power-law continuum of QSFit
import QSFit.QSORecipes.add_qso_continuum!
function add_qso_continuum!(recipe::CRecipe{T}, fp::GModelFit.FitProblem, ith::Int) where T <: WP9Type1IR
    @track_recipe
    λ = coords(getdomain(fp, ith))

    comp = QSFit.powerlaw(3000)
    #comp = QSFit.sbpl(3000)
    comp.x0.val = median(λ)
    comp.norm.val = median(values(getdata(fp, ith))) # Can't use Dierckx.Spline1D since it may fail when data is segmented (non-good channels)
    (comp.norm.val < 0)  &&  (comp.norm.val = mad(values(getdata(fp, ith))))
    comp.norm.low = comp.norm.val / 1000.  # ensure contiuum remains positive (needed to estimate EWs)
    comp.alpha.val  = -1.5 #Relevant for power law and cut-off power law
    comp.alpha.low  = -5 #Relevant for power law and cut-off power law
    comp.alpha.high =  5 #Relevant for power law and cut-off power law
    #comp.delta.val = 0.001 #Relevant for SBPL
    #comp.delta.fixed = true #Relevant for SBPL
    getmodel(fp, ith)[:QSOcont] = comp
    push!(getmodel(fp, ith)[:Continuum].list, :QSOcont)
end

function init_recipe!(recipe::CRecipe{T}) where T <: WP9Type1IR
    @track_recipe
    @invoke init_recipe!(recipe::CRecipe{<: Type1})
    recipe.line_component = QSFit.SpecLineGauss # type of line profile to use
    recipe.wavelength_range=[1216, 20000]
    recipe.host_template_range = [4000, 30000] #wavelength range to consider for the host templates to be fit
    # recipe.host_template[:ref_wavelength] = 8000.
    recipe.min_spectral_coverage[:Ironuv]= 0.1
    recipe.min_spectral_coverage[:Ironopt]= 0.1
end


function lines_dict(recipe::CRecipe{T}) where T <: WP9Type1IR
    @track_recipe

    out = QSFit.SpecLineSet()
    add_line!(recipe, out, :Lya         , QSFit.BroadLine)
    add_line!(recipe, out, :NV_1241     , QSFit.NarrowLine)
    add_line!(recipe, out, :OI_1306     , QSFit.BroadLine)
    add_line!(recipe, out, :CII_1335    , QSFit.BroadLine)
    add_line!(recipe, out, :SiIV_1400   , QSFit.BroadLine)
    add_line!(recipe, out, :CIV_1549    , QSFit.BroadLine)
    add_line!(recipe, out, :HeII_1640   , QSFit.BroadLine)
    add_line!(recipe, out, :OIII_1664   , QSFit.BroadLine)
    add_line!(recipe, out, :AlIII_1858  , QSFit.BroadLine)
    add_line!(recipe, out, :CIII_1909   , QSFit.BroadLine)
    add_line!(recipe, out, :CII_2326    , QSFit.BroadLine)
    # add_line!(recipe, out, :l2420p0     , QSFit.BroadLine)
    add_line!(recipe, out, :MgII_2798   , QSFit.BroadLine)
    add_line!(recipe, out, :NeV_3345)
    add_line!(recipe, out, :NeV_3426)
    add_line!(recipe, out, :OII_3727)
    add_line!(recipe, out, :NeIII_3869)
    add_line!(recipe, out, :Hd          , QSFit.BroadLine)
    add_line!(recipe, out, :Hg          , QSFit.BroadLine)
    add_line!(recipe, out, :OIII_4363)
    add_line!(recipe, out, :HeII_4686   , QSFit.BroadLine)
    add_line!(recipe, out, :Hb          , QSFit.BroadLine)
    add_line!(recipe, out, :OIII_4959)
    add_line!(recipe, out, :OIII_5007)
    add_line!(recipe, out, :HeI_5876    , QSFit.BroadLine)
    add_line!(recipe, out, :OI_6300)
    add_line!(recipe, out, :OI_6364)
    #add_line!(recipe, out, :NII_6549)
    add_line!(recipe, out, :Ha          , QSFit.BroadLine, QSFit.VeryBroadLine)
    # add_line!(recipe, out, :NII_6583)
    add_line!(recipe, out, :SII_6716)
    add_line!(recipe, out, :SII_6731)
    
    add_line!(recipe, out, :Paa         , QSFit.BroadLine)
    add_line!(recipe, out, :Pab         , QSFit.BroadLine)
    add_line!(recipe, out, :Pag         , QSFit.BroadLine)
    add_line!(recipe, out, :Pad         , QSFit.BroadLine)
    add_line!(recipe, out, :OI_8448     , QSFit.BroadLine)
    add_line!(recipe, out, :SIII_9532)
    add_line!(recipe, out, :HeI_10832   , QSFit.BroadLine)
    # add_line!(recipe, out, :SIII_9069)
    # add_line!(recipe, out, :OI_11289, NarrowLine, BroadLine)
    # add_line!(recipe, out, :OI_13168, NarrowLine)
    # add_line!(recipe, out, :CaII_8500, NarrowLine)
    # add_line!(recipe, out, :CaII_8544, NarrowLine)
    # add_line!(recipe, out, :CaII_8664, NarrowLine)
    # add_line!(recipe, out, :HI_9231, NarrowLine, BroadLine)
    # add_line!(recipe, out, :HI_10052, NarrowLine, BroadLine)
    # add_line!(recipe, out, :HI_10941, NarrowLine, BroadLine)
    # add_line!(recipe, out, :HI_12821, NarrowLine, BroadLine)
    # add_line!(recipe, out, :HI_18756, NarrowLine, BroadLine)
    # add_line!(recipe, out, :HI_19450, NarrowLine, BroadLine)
    # add_line!(recipe, out, :SVIII_9913, NarrowLine)
    # add_line!(recipe, out, :SIX_12523, NarrowLine)
    # add_line!(recipe, out, :CI_9826, NarrowLine)
    # add_line!(recipe, out, :CI_9852, NarrowLine)
    # add_line!(recipe, out, :HeII_10126, NarrowLine, BroadLine)
    # add_line!(recipe, out, :HeII_11629, NarrowLine)
    # add_line!(recipe, out, :SII_10289, NarrowLine)
    # add_line!(recipe, out, :SII_10323, NarrowLine)
    # add_line!(recipe, out, :NI_10400, NarrowLine)
    # add_line!(recipe, out, :NI_10410, NarrowLine)
    # add_line!(recipe, out, :PII_11471, NarrowLine)
    # add_line!(recipe, out, :PII_11886, NarrowLine)
    # add_line!(recipe, out, :SiX_14304, NarrowLine)
    # add_line!(recipe, out, :SiVI_19650, NarrowLine)
    
    #haskey(out, :NII_6549)  &&  delete!(out, :NII_6549)
    #haskey(out, :NII_6583)  &&  delete!(out, :NII_6583)
    #haskey(out, :OIII_5007_bw)  &&  delete!(out, :OIII_5007_bw)
    #delete!(out, :Ha_na)
    #delete!(out, :Ha_bb)
    #delete!(out, :Hb_na)
    #delete!(out, :MgII_2798_na)
    
    #add_line!(recipe, out, :FeII_9229, NarrowLine)
    #add_line!(recipe, out, :FeII_12570, NarrowLine)
    #add_line!(recipe, out, :FeII_12791, NarrowLine)
    #add_line!(recipe, out, :FeII_13209, NarrowLine)
    #add_line!(recipe, out, :FeII_15338, NarrowLine)
    #add_line!(recipe, out, :FeII_16439, NarrowLine)
    #add_line!(recipe, out, :FeII_16773, NarrowLine)
    #add_line!(recipe, out, :FeII_8929, NarrowLine)
    #add_line!(recipe, out, :FeII_9125, NarrowLine)
    #add_line!(recipe, out, :FeII_9134, NarrowLine)
    #add_line!(recipe, out, :FeII_9178, NarrowLine)
    #add_line!(recipe, out, :FeII_9180, NarrowLine)
    #add_line!(recipe, out, :FeII_9199, NarrowLine)
    return out
end

import QSFit.QSORecipes.add_iron_opt!
function add_iron_opt!(recipe::CRecipe{<: WP9Type1IR}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe    
    @invoke add_iron_opt!(recipe::CRecipe{<: Type1}, fp, ith)
    model = QSFit.QSORecipes.getmodel(fp, ith)
    if haskey(model, :Ironoptna)
        model[:Ironoptna].norm.val = 0.
        model[:Ironoptna].norm.fixed = true
    end
end

import QSFit.QSORecipes: add_patch_functs!
function add_patch_functs!(recipe::CRecipe{<: WP9Type1IR}, fp::GModelFit.FitProblem, ith::Int)
    @track_recipe
    model = getmodel(fp, ith)
    if  haskey(model, :Hb_br)  &&
        haskey(model, :Hb_na)
        model[:Hb_na].fwhm.high = 0.3
        model[:Hb_na].fwhm.low  = 0.
        model[:Hb_na].fwhm.val  = 0.1
        model[:Hb_na].fwhm.patch = @fd (m, v) -> v * m[:Hb_br].fwhm
    end

    if  haskey(model, :Ha_br)  &&
        haskey(model, :Ha_na)
        model[:Ha_na].fwhm.high = 0.3
        model[:Ha_na].fwhm.low  = 1. / 20
        model[:Ha_na].fwhm.val  = 0.1
        model[:Ha_na].fwhm.patch = @fd (m, v) -> v * m[:Ha_br].fwhm
    end
end



end # module LinesInTheSky
