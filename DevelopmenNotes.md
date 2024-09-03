# Development notes

This package in the current form is under heavy restructuring.
Here are notes on not-obvious facts.

## The Problem
Most of the logic in the original PMFRG package
does not need OrdinaryDiffEq to be tested,
but nonetheless precompilation 
(and test of the logic that includes solve)
need it to be precompiled.

This means that every time one does 
```
] test PMFRG
```

Tests might take a lot of time because PMFRG
and the dependencies need to be recompiled,
and OrdinaryDiffEq is notoriously long to compile.

OrdinaryDiffEq was recently split into many different packages
that contain different methods,
so that it is possible to include only the sub-package 
that contains the method of interest,
and this saves some compilation time.
However this approach seem to save only 3m17s 
out of 18m52s, 
(on GitHub Actions)
so it is not a perfect solution.

## Proposed solution

Split the PMFRG.jl package into PMFRGCore.jl,
that does not need anything 
of the Differential Equation ecosystem 
to be tested,
that can be tested quickly 
and will contain the vast majority of the code,
and another package PMFRGSolve.jl,
which contains the minimum necessary
that is needed to couple PMFRGCore.jl
with the differential equation solvers
(thingsin OrdinaryDiffEq subpackages,
and solver Callbacks).

Then create an umbrella package
PMFRG.jl that uses both and offers 
the same interface as the original PMFRG.jl.

Note: since the two sub packages need to be in sync,
they need to belong to the same git repository.

For convenience, we have all 3 in the same repo.

**Point for discussion**: can't we just have PMFRGSolve.jl == PMFRG.jl?


## Dependencies

Since we want to iterate on the subpackages 
PMFRGCore.jl and PMFRGSolve.jl,
and they are not in any index,
they need to be added as dev dependencies 
to the umbrella package PMFRG.jl.

Also, PMFRGSolve.jl might need PMFRGCore.jl,
as a dev package.

so I have done
```
] activate PMFRGSolve.jl
] dev ./PMFRGCore.jl/
] activate . # meaning, the umbrella package
] dev ./PMFRGCore.jl/ ./PMFRGSolve.jl/ 
```

**CRUCIAL**: add them both with the same dev command,
otherwise it will fail saying that one of the two is not in an index.
This is likely related to the fact that PMFRGSolve.jl
itself depends already on PMFRGCore.jl. 
Maybe this is related to this **Point for discussion**: 
can't we just have PMFRGSolve.jl == PMFRG.jl?
