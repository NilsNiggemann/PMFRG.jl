# Notes for the MPI parallelization

## Profiling

The code has been profiled
using the first example in the README,
which has been wrapped in a script 
(`benchmarks/profile.jl`)

### Multithreaded Benchmarking

#### Laptop
Benchmarking run on my own laptop (4-Core i5-1145G7 @ 2.60 GHz).
Some other programs were idling when I was running the benchmarks.
I used a modified version 
of the  `Example/generating_data.jl` script
as a benchmark (see `generating_data_benchmark.jl` in this directory)
and noted down the time measurement
that was printed out by the `@time` macro inside `launchPMFRG!` in `Solver.jl`.

Using these parameters:
```

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=30, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)

```

With multiple threads we have:

+----------+----------+
| Nthreads | Time (s) |
+----------+----------+
|        1 |     21.3 |
|        2 |     14.6 |
|        4 |      9.4 |
|        8 |      7.3 |
|       16 |     10.7 |
+----------+----------+


Switching to `N=50`:

```
Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=50, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)
```

With multiple threads we have 

+----------+----------+
| Nthreads | Time (s) |
+----------+----------+
|        1 |    142.7 |
|        2 |     91.3 |
|        4 |     66.0 |
|        8 |     60.9 |
|       16 |     57.3 |
+----------+----------+

