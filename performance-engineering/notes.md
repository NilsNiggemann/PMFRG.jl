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

#### Scaling with N (76 threads, NLen=14)

Wondering 
Take 1 (cpuonly partition)

| N  | Time (s)| Retake 1 |
|----|---------|----------|
| 20 |    10.5 |     14.6 |
| 21 |    12.4 |     15.9 | 
| 22 |    14.7 |     20.8 |
| 23 |    17.1 |     23.4 |
| 24 |    20.0 |     27.5 |
| 25 |    23.1 |     34.3 |

Recompiling and relaunching (accelerated partition)

| N  | Time (s)|
|----|---------|
| 20 |    14.8 |
| 21 |    15.9 |
| 22 |    20.7 |
| 23 |    23.3 |
| 24 |    27.2 |
| 25 |    33.8 |

## MPI-ization roadmap

1. Parallelization of getXBubble!
  1a. MPI-ize getXBubble! 
     1. Create regression tests using Recorder.jl.
        Make sure that the recorded test cases are significant
        (e.g., the recorded arrays are not full of zeros, for example) [DONE]
     2. Create getXBubblePartition! that works only on a part of the array,
        use that inside getXBubble!, verify correctness 
        with the recorded test cases.
        getXBubblePartition should contain all the possible functionality
        that does not require MPI calls. [DONE]
        OPEN ISSUE: how much logic should be outside 
        of getXBubblePartition! ?
        This is related to 1b.2.
     3. map between partitions of is,it,iu ranges and mpi ranks. [DONE]
        (computing ranges might be a little expensive, 
        if so better to cache them,
        perhaps with Memoize.jl? Measure first)
        
  1b. Add MPI setup to original program
    1. Create functions to broadcast portions of arrays
       (already using MPI). 
       This is trivial but a working example must be produced [DONE]
    2. Add MPI initialization and finalization in test script
       using only non MPI-functions.
       Conflicts on external resources must be solved somehow
       (e.g., making so that each process writes in his own directory,
       named after the process itself, or some other solution).
       No MPI-ized functions should be implemented [DONE on one example]
    3. Implement MPI version of getXBubble!, called getXBubbleMPI!,
       as a drop-in replacement for getXBubble!,
       that should work without any changes 
       in the script with the MPI setup. [DONE?]
       1. Create mpi test for getXBubbleMPI!,
          using regression test case data from getXBubble! [DONE]
       2. MPI.Initialized() can be called 
          to determine which version
          (either getXBubbleMPI or getXBubble)
          can be called. 
     
  1c. CRUCIAL: Add all the tests to the test suite
      (a bit tricky with the mpi ones) [DONE]

  2. Using MPI.Initialized(), 
     move the logic of getXBubbleMPI! inside getXBubble!
     in an `if MPI.Initialized()` conditional
     and remove getXBubbleMPI!. [DONE]

  3. Allow superposition of communication and computation
     by assigning more portions to each rank (review 1a.3)
     and either use one `MPI_Ibcast` per portion 
     or one `MPI_Isend` / `MPI_Irecv` per portion for rank.
  
