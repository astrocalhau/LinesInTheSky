# WP9Type1IRRecipe
This is hopefully something that will eventually become the recipe for type I AGN to use with QSFIT for Euclid's WP9 group.

At present, it is only a copy of Giorgio Calderone's Type2IR recipe, with the references to type II changed to type I and a couple of comments of mine where I try to understand what the elements do and how to build a recipe from an existing one.

Whatever comes out of this, if it works, it works because of Giorgio and Anna's help. If it doesn't work, I probably messed up somewhere.

## Usage
The steps for usage of this recipe follow closely to the original work by Giorgio Calderone.
- Install Julia (https://julialang.org/downloads/);
- Download this folder;
- Start a Julia session specifying the path to WP9Type1IRRecipe as the current environment: julia --project=path/to/WP9Type1IRRecipe;
- If this is the first time running the recipe on this session, then:
  - Download and compile the dependencies by typing the following code in the Julia terminal:
      using Pkg
      @assert Pkg.project().name == "WP9Type1IRRecipe"
      Pkg.instantiate()
- to test the package run "type1_recipe_test.jl". This can be done directly on the step for starting a julia session by writing instead: julia --project=path/to/WP9Type1IRRecipe type1_recipe_test.jl
  The test file should take care of the instantiation by itself.
- the separate steps to run the analysis are:

          using QSFit, GModelFitViewer

          using WP9Type1IRRecipe

          spec = Spectrum(Val(:ASCII), label="yourlabel", columns=[1,2],
                          "path/to/input/spectra.txt")

          recipe = CRecipe{WP9Type1IR}(redshift=0.355 (replace as needed))

          res = analyze(recipe, spec)

          viewer(res)

          GModelFitViewer.serialize_html(filename="path/to/output.html", res, keep=true)

## Beta version, I guess.

