# Notes for the MPI parallelization

## Profiling

The code has been profiled
using the first example in the README,
which has been wrapped in a script 
(`benchmarks/profile.jl`)

### Multithreaded Benchmarking

#### Toy problem with NLen=5
##### Laptop
Benchmarking run on my own laptop (4-Core i5-1145G7 @ 2.60 GHz).
Some other programs were idling when I was running the benchmarks.
I used a modified version 
of the  `Example/generating_data.jl` script
as a benchmark (see `generating_data_benchmark.jl` in this directory)
and noted down the time measurement
that was printed out by the `@time` macro inside `launchPMFRG!` in `Solver.jl`.

Using these parameters:
```
System = getSquareLattice(5, [1,0.1])

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

| Nthreads | Time (s) |
|----------|----------|
|        1 |     21.3 |
|        2 |     14.6 |
|        4 |      9.4 |
|        8 |      7.3 |
|       16 |     10.7 |


Switching to `N=50`:

```
System = getSquareLattice(5, [1,0.1])

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

| Nthreads | Time (s) |
|----------|----------|
|        1 |    142.7 |
|        2 |     91.3 |
|        4 |     66.0 |
|        8 |     60.9 |
|       16 |     57.3 |

##### HoreKa 
NOTE: This case was never meant to run on a supercomputer.
Benchmarking on a single node HoreKa, 76-core machine
(Intel Xeon Platinum 8368 @ 2.40GHz).
Not using any pinning to any core.

```
System = getSquareLattice(5, [1,0.1])

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

| Nthreads | Time (s) |
|----------|----------|
|        1 |    178.9 |
|        2 |    100.6 |
|        4 |     58.9 |
|        8 |     38.7 |
|       19 |     21.5 |
|       38 |     14.5 |
|       76 |     11.9 |
|      152 |     17.2 |


#### "Real" model with NLen=14

This combination of parameters is possibly more realistic
and should give a better idea of the performance of the real code.

```julia
NLen = 14 
J1 = 1
J2 = 0.1
couplings = [J1, J2] 

System = getSquareLattice(NLen, couplings) 
Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N=25, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy=1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
    MinimalOutput=true,
)
```

| Nthreads | Time - floops (s) | Time - threads (s) |
|----------|-------------------|--------------------|
|        1 |             408.2 |              408.0 |
|        2 |             216.8 |              208.4 |
|        3 |                   |              141.7 |
|        4 |             122.3 |              108.7 |
|        5 |                   |               88.0 |
|        6 |                   |               73.9 |
|        7 |                   |               64.4 |
|        8 |              73.4 |               57.3 |
|       10 |                   |               47.4 |
|       15 |                   |               33.5 |
|       19 |              39.1 |               28.2 |
|       25 |                   |               22.6 |
|       30 |                   |               20.5 |
|       38 |              22.5 |               16.5 |
|       50 |                   |               15.1 |
|       60 |                   |               15.9 |
|       76 |              19.9 |               14.3 |
|      100 |                   |               12.7 |
|      120 |                   |               19.2 |
|      150 |                   |               27.0 |
|      152 |              26.6 |               28.7 |



It could be that this problem 
(with the current implementation) 
is too small for the node.
#### "Real" model with NLen=20 and N=50

This model is too large for profiling:

| Nthreads | Time (s) |
|----------|----------|
|      152 |    576.8 | 
|       76 |    665.7 |




