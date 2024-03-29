using PMFRG, MPI, TimerOutputs
include("mpi/MPI_Detail.jl")
import .MPI_Detail

function PMFRG.getXBubble!(Workspace, Lam, ::PMFRG.UseMPI)
    (; X, Par) = Workspace
    (; N, np_vec) = Par.NumericalParams

    if MPI.Initialized()
        nranks = MPI.Comm_size(MPI.COMM_WORLD)
        rank = MPI.Comm_rank(MPI.COMM_WORLD)

        # if (ns+nt+nu)%2 == 0	# skip unphysical bosonic frequency combinations
        #
        if (3 * np_vec[1]) % 2 == 0
            # then 1,1,1 , with parity 1, is not right
            parity = 1
        else
            parity = 0
        end

        @timeit_debug "get_ranges" all_ranges =
            MPI_Detail.get_all_ranges_stu(N, nranks, parity)
        iurange_full = 1:N
        isrange, itrange, _ = all_ranges[rank+1]
        @timeit_debug "partition" PMFRG.getXBubblePartition!(
            X,
            Workspace.State,
            Workspace.Deriv,
            Par,
            Workspace.Buffer,
            Lam,
            isrange,
            itrange,
            iurange_full,
        )


        @timeit_debug "communication" for root = 0:(nranks-1)
            isrange, itrange, iurange_restrict = all_ranges[root+1]
            iurange_abc = Par.Options.usesymmetry ? iurange_restrict : iurange_full

            MPI.Bcast!((@view X.a[:, isrange, itrange, iurange_abc]), root, MPI.COMM_WORLD)
            MPI.Bcast!((@view X.b[:, isrange, itrange, iurange_abc]), root, MPI.COMM_WORLD)
            MPI.Bcast!((@view X.c[:, isrange, itrange, iurange_abc]), root, MPI.COMM_WORLD)

            MPI.Bcast!(
                (@view X.Ta[:, isrange, itrange, iurange_full]),
                root,
                MPI.COMM_WORLD,
            )
            MPI.Bcast!(
                (@view X.Tb[:, isrange, itrange, iurange_full]),
                root,
                MPI.COMM_WORLD,
            )
            MPI.Bcast!(
                (@view X.Tc[:, isrange, itrange, iurange_full]),
                root,
                MPI.COMM_WORLD,
            )
            MPI.Bcast!(
                (@view X.Td[:, isrange, itrange, iurange_full]),
                root,
                MPI.COMM_WORLD,
            )
        end
    else
        @warn "MPI package used but not initialized" maxlog = 1
        PMFRG.getXBubblePartition!(X, State, Deriv, Par, Buffer, Lam, 1:N, 1:N, 1:N)
    end
end
