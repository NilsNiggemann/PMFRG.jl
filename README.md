# PMFRG.jl: 

`PMFRG.jl` (**P**seudo-**M**ajorana **F**unctional **R**enormalization **G**roup) is a Julia package used to compute observables for spin- $1/2$ Heisenberg models of the form

```math
H = \sum_{ij} J_{ij} \vec{S}_i \cdot \vec{S}_j
```
## Installation
`PMFRG.jl` is not in the official registry (yet?) and thus best installed via installing the private registry "JuliaPMFRG":
```
(@v1.9) pkg> registry add https://github.com/NilsNiggemann/JuliaPMFRGRegistry.git
```
This only needs to be done once on every machine. Then, `PMFRG.jl`, its dependencies and other helper packages for evaluation can be installed conveniently via
```
(v1.9) pkg> add PMFRG
```
Note, that is advised to create a reproducible environment for each project, see also https://pkgdocs.julialang.org/v1/environments/ .

<details>
  <summary>Installation without local registry</summary>
  
If you do not want to use the private registry, you can also install the dependencies manually.

```
(@v1.9) pkg> activate TestProject
Activating new project at `~/TestProject`

(TestProject) pkg> add https://github.com/NilsNiggemann/SpinFRGLattices.jl.git
(TestProject) pkg> add https://github.com/NilsNiggemann/PMFRG.jl.git
```
Note that in this case the dependency `SpinFRGLattices` needs to be added first for the package manager to resolve the correct version. 

</details>

## Usage
After the package is installed to a local environment it can be loaded in a Julia session with  `using PMFRG`. The following contains a minimal working example:
```
using SpinFRGLattices, PMFRG
using SpinFRGLattices.SquareLattice

NLen = 5 # Number of nearest neighbor bonds up to which correlations are treated in the lattice. For NLen = 5, all correlations C_{ij} are zero if sites i and j are separated by more than 5 nearest neighbor bonds.
J1 = 1
J2 = 0.5
couplings = [J1,J2] # Construct a vector of couplings: nearest neighbor coupling is J1 (J2) and further couplings to zero. For finite further couplings simply provide a longer array, i.e [J1,J2,J3,...]

System = getSquareLattice(NLen,couplings) # create a structure that contains all information about the geometry of the problem. 

Par = Params( #create a group of all parameters to pass them to the FRG Solver
    System, # geometry, this is always required
    OneLoop(), # method. OneLoop() is the default
    T=0.5, # Temperature for the simulation.
    N = 10, # Number of positive Matsubara frequencies for the four-point vertex.
    accuracy = 1e-3, #absolute and relative tolerance of the ODE solver.
    # For further optional arguments, see documentation of 'NumericalParams'
)

mainFile = "temp/"*PMFRG.generateFileName(Par,"_testFile") # specify a file name for main Output
flowpath = "temp/flows/" # specify path for vertex checkpoints

Solution,saved_values = SolveFRG(Par,MainFile = mainFile ,CheckpointDirectory = flowpath,method = DP5(),VertexCheckpoints = [],CheckPointSteps = 3 )
```
For further options, the documentation of `SolveFRG` or `NumericalParams` can be helpful. Note that if no `MainFile` is specified, then no output is written.
### Beyond oneloop FRG
If no further specification is made, the standard oneloop implementation will be used. In most cases, the oneloop approximation should be appropriate. The Twoloop implementation can be used via specifying the parameters with the additional argument `Params(System,TwoLoop();kwargs...)`.
## Usage on SLURM Clusters
To start a run on a cluster, we need both a batch script which requests the required resources and a job-script that runs the calculation. An example is found in the file ´Example/Slurm_example.jl´ in this repository.
Since Julia ignores everything after `#`, we can also place the batch script and the job script in the same file. To run a job array with indices `1-10`, we only need to run `sbatch --array=1-10 Example/Slurm_example.jl` in the cluster's terminal on the login node. The indices of the job array are passed as arguments, so they can be accessed via `i_arg = parse(Int, ARGS[1])`. 
One can also sweep over several parameters at once by linearly indexing an array that contains all the combinations of parameters:
```
Trange = 0.2:0.1:2.
J2range = 0.025:0.025:0.175

jobsarray = [(t,j) for t in Trange,j in J2range]
T,J2 = jobsarray[i_arg]
```
## Fine-tuning lattice couplings
We might also want to look at lattices where the couplings are not based on distance. The couplings between each **symmetry inequivalent** pair of sites can also be set individually for example
```
setCoupling!(System,1,Rvec(1,0,2),0.3)
```
sets the coupling from reference site 1 (`Rvec(0,0,1)`) to another site located at Rvec(1,0,2) to a value of `0.3`. 
Generally, `Rvec(n1,n2,nb)` means the site is located at 
 $` n_1 \vec{a}_1 +n_2\vec{a}_2+ \vec{b}_{nb}`$ where $`\vec{a}_i`$ are lattice vectors and $`\vec{b}_{nb}`$ is the $`n_b`$ 'th basis vector.
Analogously a site in a generic 3D lattice is given by `Rvec(n1,n2,n3,nb)`.
**attention:** Setting the couplings can **never** break any of the hard-coded lattice symmetries. The function `setCoupling!` will search for the corresponding symmetry inequivalent pair instead. For the square lattice, which is here implemented using all its mirror symmetries the following two lines lead to equivalent results:

```
setCoupling!(System,1,Rvec(-1,0,1),0.3)
setCoupling!(System,1,Rvec(1,0,1),0.3)
```
If you need to set couplings that do not correspond to a lattices symmetry, the only option is to implement a new version of the lattice with reduced symmetry, which will increase the numerical complexity.
If you are implementing a new lattice, or experimenting with complicated couplings, it can be helpful to use the package `FRGLatticePlotting`:

```
using SpinFRGLattices, SpinFRGLattices.SquareLattice, FRGLatticePlotting
System = getSquareLattice(5,[1.,0])
setCoupling!(System,1,Rvec(-1,0,1),0.3)
testGeometry(System)
plotSystem(System,Basis)
```
This will perform a bunch of standardized unit tests that should be fulfilled by a correct implementation. The last part will display a plot of the lattice indicating symmetry-inequivalent pairs as well as their couplings.  
Another helpful trick is to plot a Fourier transform of the couplings (the list of couplings takes the same form as the list of susceptibilities, see further below). If a lattice is implemented correctly this plot has to obey all symmetries that you would also expect for the susceptibility.
## Output
Typically the main output consists of a single HDF5 file that can be loaded using the HDF5 library.
It contains a number of observables during the course of the flow. If no group is specified in `SolveFRG` then the temperature is used for the main Group. This way, HDF5 files can be easily merged into a single file for several temperatures. Example:
```
using HDF5
julia> f = h5open(fn_PM)
🗂️ HDF5.File: (read-only) "PMFRG_SquareLattice_NLen=5_N=10_l1.h5"
├─ 📂 0.5
│  ├─ 🔢 Chi
│  ├─ 🔢 Chi_nu
│  ├─ 🔢 Lambda
│  ├─ 🔢 MaxVa
│  ├─ 🔢 MaxVb
│  ├─ 🔢 MaxVc
│  ├─ 🔢 N
│  ├─ 🔢 NLen
│  ├─ 🔢 NUnique
│  ├─ 🔢 T
│  ├─ 🔢 f_int
│  └─ 🔢 gamma
julia> chi = read(f["0.5/Chi"]);
julia> wholefileasDict = read(f);
close(f)
```
For longer runs, it is advisable to specify a `CheckPointDirectory`, where regular Checkpoints of the Solution state are stored. From each checkpoint, the Solver can be re-started. Also, the data can be used to for further anaysis, if complete vertex information should be required.
Note that a single Checkpoint may take up several Gigabytes of storage (for cluster-sized problems) so they should be saved to the scratch directory.
## Evaluation
Since susceptibilities are returned as a list according to pre-selected symmetry inequivalent pairs, the library [`SpinFRGLattices.jl`](https://github.com/NilsNiggemann/SpinFRGLattices.jl)) has to be used for evaluation. To compute Fourier transforms (and other things), it is helpful to use the package [`PMFRGEvaluation`](`https://github.com/NilsNiggemann/PMFRGEvaluation.jl`). As an example, the following code plots the magnetic susceptibility for the data generated above.
```
using HDF5, PMFRGEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course

Lattice = LatticeInfo(System,SquareLattice)

chi_ΛR = h5read(mainFile,"0.5/Chi")
chi_R = chi_ΛR[:,end]

chi = getNaiveLatticeFT(chi_R,Lattice)
# or: chi = getLatticeFFT(chi_R, Lattice) to use FFT instead
k = LinRange(-2pi,2pi,100)

chik = [chi(x,y) for x in k, y in k]

heatmap(k,k,chik)

```

### Notes
- The object `chi = getNaiveLatticeFT(...)` corresponds to the average over all sublattices $\chi(\vec{k}) = \frac{1}{N_\textrm{Cell}} \sum_{a,b} \chi_{ab}(\vec{k})$. 
- For lattices with more than one site per unit cell, one can also evaluate the correlations between sublattice 1 and 2 as `chi[1,2](kx,ky,kz)`. It is also possible to add two susceptibilities (before evaluating them at some momentum) such as `chi_new = chi[1,2] + chi[2,1]`.

- If you are going to evaluate the Fourier transform in the whole Brillouin zone, it might be more efficient to use the Fast Fourier Transform (FFT) instead. This can be done by using `chi = getLatticeFFT(chi_R, Lattice)` instead. Note that continuous frequencies are here obtained by zero-padding the data which is less accurate than the naive method. The padding default can be changed for example to include at least $128$ sites in each direction as `chi = getLatticeFFT(chi_R, Lattice, 128)`.


## More Examples
A more thorough set of examples is found in the `Examples` folder of this repository. For code reuse, the dependencies of [`PMFRGEvaluation`](`https://github.com/NilsNiggemann/PMFRGEvaluation.jl`) are split into several subdependencies. To try out the examples, activate the project environment with `]activate Example` and download all dependencies with `]instantiate`.

It is a good practice to set up a new evaluation environment for each project. If you use the same environment for everything, you might not be able to reproduce plots you made a while ago, because the plotting package or the evaluation package may have changed.

## MPI parallelization
By default, the package uses `Threads` 
to use all the cores of a single node of a HPC machine.
In order to use more than one node,
this package includes an [extension](https://pkgdocs.julialang.org/v1/creating-packages/#Conditional-loading-of-code-in-packages-(Extensions))
that can use [MPI](https://en.wikipedia.org/wiki/Message_Passing_Interface)
via the [MPI.jl wrapper](https://juliaparallel.org/MPI.jl/stable/).

In order to activate the extension the following conditions are necessary:
- the Julia version needs to be >= 1.9 (the package extension mechanism is not implemented in older versions)
- the MPI package needs to be added to the active environment 
Then, the Julia script can be launched using [`mpiexecjl`](https://juliaparallel.org/MPI.jl/stable/usage/#Julia-wrapper-for-mpiexec).

An example of usage is in `ext/PMFRGMPIExt/test/MPITest/generate_data_example.mpi.jl`
(only the Julia code is present:
ideas for a wrapping slurm/shell script 
can be taken from `Example/Slurm_example.jl`).
Some notes:
- It is necessary to initialize MPI at the beginning by calling `MPI.Init()`,
  (and call `MPI.Finalize()` at the end).
  If `MPI.Init()` is not called, the MPI processes will not be able to communicate.
  The code will launch a warning in this case.
- In order for the MPI version of the code to be called, 
  please pass a `UseMPI()` singleton argument to the `SolveFRG` function
  (see example below).
- Any MPI launcher used will launch the same Julia script on different processes.
  For this reason, it is important to make sure 
  that either the output of the processes 
  (in MPI terminology, *ranks*)
  is written on different files (which should be identical),
  or only one MPI process writes the outputs to file.
  In the following code snippet,
  only the master rank will write the output to disk:
  
  
``` julia
using MPI

MPI.Init()
rank = MPI.Comm_rank(MPI.COMM_WORLD)

[...]

if rank == 0
    # specify a file name for main Output
    mainFile = PMFRG.generateFileName(Par, "_testFile")
    # specify path for vertex checkpoints
    flowpath = "flows/"
else
    # disable file output for other ranks
    mainFile = nothing
    flowpath = nothing
end


Solution, saved_values = SolveFRG(
    Par,
    UseMPI(),
    MainFile = mainFile,
    CheckpointDirectory = flowpath,
    method = DP5(),
    VertexCheckpoints = [],
    CheckPointSteps = 3,
)


MPI.Finalize()
```

   The output from all the ranks (e.g, the files in all the directories) should match.

## Implementing your own lattices
Of course you will eventually have to implement lattices which are not included already, change the couplings, or even remove symmetries. As long as you feed a valid geometry struct from [`SpinFRGLattices.jl`](https://github.com/NilsNiggemann/SpinFRGLattices.jl) to the FRG code (which is quite minimalistic), it should not be necessary to make direct changes to the library. [`SpinFRGLattices.jl`](https://github.com/NilsNiggemann/SpinFRGLattices.jl)) contains mostly helper functions to make your life easier. Documentation of how to use it to implement new lattices is found soon in the repository.
## See also

- Lattice implementations: [`SpinFRGLattices.jl`](https://github.com/NilsNiggemann/SpinFRGLattices.jl)
- Evaluation: [`PMFRGEvaluation.jl`](https://github.com/NilsNiggemann/PMFRGEvaluation.jl)
- visualization of lattices: [`FRGLatticePlotting.jl`](https://github.com/NilsNiggemann/FRGLatticePlotting.jl)
### PMFRG methodological publications
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.103.104431
- https://scipost.org/SciPostPhys.12.5.156
### Other packages of interest
- For zero temperature calculations: [`PFFRGSolver.jl`](https://github.com/dominikkiese/PFFRGSolver.jl)
