# PMFRG.jl
`PMFRG` is a Julia package used to compute observables for spin-$1/2$ Heisenberg models of the form
$$
H = \sum_{ij} J_{ij} \vec{S}_i \cdot \vec{S}_j
$$
## Installation
 Currently, you need to have access to this repository to install the package. If you can read this then you have access and the following should work. It is advise to create a reproducible environment for each project, see also https://pkgdocs.julialang.org/v1/environments/ . Note that to install `PMFRG.jl` the package `SpinFRGLattices.jl` is also required. To install both, just paste the following into the Julia REPL:
```
]add git@gitlabph.physik.fu-berlin.de:niggeni/spinfrglattices.jl.git, "git@gitlabph.physik.fu-berlin.de:niggeni/pmfrg.jl.git"
```
If ssh authentication is not possible, the package can also be installed using https instead. The following will require a password to the repository:
```
]add https://gitlabph.physik.fu-berlin.de/niggeni/spinfrglattices.jl.git, https://gitlabph.physik.fu-berlin.de/niggeni/pmfrg.jl.git
```

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

## Evaluation
Since susceptibilities are returned as a list according to pre-selected symmetry inequivalent pairs, the library `SpinFRGLattices` has to be used for evaluation. To compute Fourier transforms (and other things), it is helpful to use the package `PMFRGEvaluation` (`git@gitlabph.physik.fu-berlin.de:niggeni/PMFRGEvaluation.git`). As an example, the following code plots the magnetic susceptibility for the data generated above.
```
using HDF5, PMFRGEvaluation
using CairoMakie #for plotting. You can use whatever plotting package you like of course

Lattice = LatticeInfo(System,SquareLattice)

chi_ΛR = h5read(mainFile,"0.5/Chi")
chi_R = chi_ΛR[:,end]

chi = getFourier(chi_R,Lattice)

k = LinRange(-2pi,2pi,100)

chik = [chi(x,y) for x in k, y in k]

heatmap(k,k,chik)

```
## More Examples
A more thorough set of examples is found in the Examples folder of this repository. For code reuse, the dependencies of `PMFRGEvaluation` are split into several subdependencies. To try out the examples, activate the project environment with `]activate Example` and download all dependencies with `]instantiate`.
To install packages in a new environment, with up-to-date packages, we need to manually download the private sub-repositories (I don't know, why). The example project was initialized with:

```
add git@gitlabph.physik.fu-berlin.de:niggeni/spinfrglattices.jl.git, git@gitlabph.physik.fu-berlin.de:niggeni/HDF5Helpers.jl.git,git@gitlabph.physik.fu-berlin.de:niggeni/pmfrg.jl.git, git@gitlabph.physik.fu-berlin.de:niggeni/FRGLatticeEvaluation.jl.git,git@gitlabph.physik.fu-berlin.de:niggeni/PMFRGEvaluation.git, CairoMakie
```
Probably, this can be made more convenient, possibly by switching to public github repositories (which is the plan anyway).

I recommend setting up a new evaluation environment for each project. If you use the same for everything, you might not be able to reproduce plots you made a while ago, because the plotting package or the evaluation package may have changed.