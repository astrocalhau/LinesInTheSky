# Lines in the sky

Paper: (TODO: add arXiv, ADS link)

## Usage
- Install Julia (https://julialang.org/downloads/);
- Download and unzip the [code](https://github.com/astrocalhau/LinesInTheSky/archive/refs/heads/main.zip);

- Download and unzip data (TODO: add link) in the `LinesInTheSky` folder;

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

- Run the analisys with: `include("run.jl")`.

- The output catalog will be written in `results/QSFIT_RESULTS.fits` and individual output files will be located in `results/HTML`.
