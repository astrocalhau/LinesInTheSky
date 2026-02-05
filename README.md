# Lines in the sky

Paper: (TODO: add arXiv, ADS link)

## Usage
- Install Julia (https://julialang.org/downloads/);
- Download and unzip the [code](https://github.com/astrocalhau/LinesInTheSky/archive/refs/heads/main.zip);

- Download the [data](https://drive.google.com/file/d/1p-4CfG_O5aWJwiIBUAkyff2sg8mwAzRF/view?usp=drive_link) in the `LinesInTheSky` folder, then untar using the command:
```
tar xvf input.tar.gz
```

- In the `LinesInTheSky` folder start julia with the following command line:
```
julia --project=. -t auto
```

- Install code depencies with:
```
using Pkg
Pkg.instantiate()
```
(note: this step is necessary only the first time you run the code)

- Depending on whether you want to analyze Euclid or DESI spectra, run the analisys with:
```
include("run_Euclid.jl")
input_path, output_path, catalog, results = run_Euclid()
```
or
```
include("run_DESI.jl")
input_path, output_path, catalog, results = run_DESI()
```

- The output catalog will be written in `results_<instrument>/QSFIT_RESULTS.fits` and individual output files will be located in `results_<instrument>/JSON`, where `instrument` is `Euclid` or `DESI` according to the choice made above.

- To reproduce most of the plots in the paper, run the following command in the bash terminal. You must be in the LinesInTheSky folder.
```
julia --project=. -t auto Make_plots.jl
```

More information on the running of QSFit and the production of user-defined recipes can be found in QSFit's official [page](https://gcalderone.github.io/QSFit.jl/)
