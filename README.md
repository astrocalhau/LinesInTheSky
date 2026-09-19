# Lines in the sky

Paper: 

ArXiV

ADS

If you are only interested in the final output tables, you can get them [here](https://drive.google.com/file/d/1txXBGhrfv7xWRckxJOFjso4857_3k6Na/view?usp=sharing)

## Usage
- Install Julia (https://julialang.org/downloads/);
- Download and unzip the [code](https://github.com/astrocalhau/LinesInTheSky/archive/refs/heads/main.zip);

- Download the [data](https://drive.google.com/file/d/1zl_P44FxsKwGNf14FYkHrglCg1jNl6CL/view?usp=sharing) in the `LinesInTheSky` folder, then untar using the command:
```
tar xvf input.tar.gz
```

- The first time you run the code:
	- In the `LinesInTheSky` folder start julia with the following command line:
	```
	julia --project=.
	```

	- Install code depencies with:

	```
	using Pkg
	Pkg.instantiate()
	```

	- Install the required dust maps with:
	```
	using DustExtinction
	dustmap = SFD98Map()
	```
	You will receive a prompt to install the dust maps, which you should accept. This only happens the first time the package is run.
	
	- Exit julia with:
	```
	exit()
	```

- In the LinesInTheSky folder start julia with multithreading with the following command line (This is also the starting point for all subsequent runs of the code):
```
julia --project=. -t auto
```

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
