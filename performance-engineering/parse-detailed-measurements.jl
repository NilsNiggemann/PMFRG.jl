
import Statistics: mean, std
import DataFrames: DataFrame, select, filter, innerjoin
using Plots


function parse_file(filename)
    open(filename) do file
        lines = readlines(file)

        total_time =
            [l for l in lines if contains(l, "seconds")][end-1] |>
            split |>
            first |>
            x -> parse(Float32, x)

        detailed_measurement_lines = [l for l in lines if contains(l, "tag:64")]

        timed_func_names = Set([split(l)[1] for l in detailed_measurement_lines])

        measurements = [
            (
                func_name = func_name,
                measurements = [
                    parse(Float32, split(l)[3]) for
                    l in detailed_measurement_lines if split(l)[1] == func_name
                ],
            ) for func_name in timed_func_names
        ]


        measurement_stats = [
            (
                funcname = row.func_name,
                total = sum(row.measurements),
                total_std = std(row.measurements) / mean(row.measurements) *
                            sum(row.measurements),
                mean = mean(row.measurements),
                std = std(row.measurements),
            ) for row in measurements
        ]

        push!(
            measurement_stats,
            (
                funcname = "Total Time",
                total = total_time,
                total_std = 0,
                mean = total_time,
                std = 0,
            ),
        )

        push!(timed_func_names, "Total Time")

        return (stats = measurement_stats, func_names = timed_func_names)
    end
end

get_nthreads(filename) = parse(Int32, match(r"detailed-timings?-(\d+)\.out", filename)[1])

function parse_list(filenames)

    all_parses =
        [(fname = fname, stats_and_funcnames = parse_file(fname)) for fname in filenames]

    return (
        fnames = all_parses[1].stats_and_funcnames.func_names,
        rows = [
            (
                Nthreads = get_nthreads(parse.fname),
                FuncName = row.funcname,
                Total = row.total,
                TotalStd = row.total_std,
                Mean = row.mean,
                Std = row.std,
            ) for parse in all_parses for row in parse.stats_and_funcnames.stats
        ],
    )
end


function make_dfs_mean(filenames)

    (; fnames, rows) = parse_list(filenames)

    df = rows |> DataFrame |> sort

    dfs = [
        (
            funcname = funcname,
            df = filter(row -> row["FuncName"] == funcname, df) |>
                 x -> select(x, "Nthreads", "Mean", "Std"),
        ) for funcname in fnames
    ]
end

function make_dfs(filenames)

    (; fnames, rows) = parse_list(filenames)

    df = rows |> DataFrame |> sort

    dfs = [
        (
            funcname = funcname,
            df = filter(row -> row["FuncName"] == funcname, df) |>
                 x -> select(x, "Nthreads", "Total", "TotalStd"),
        ) for funcname in fnames
    ]
end


function plot_dfs(dfs, savename)

    markers =
        filter(
            m -> m ∈ Plots.supported_markers() && m ∉ [:cross, :xcross],
            Plots._shape_keys,
        ) |> Iterators.cycle

    m, s = iterate(markers)
    plot()

    for d in dfs
        p = scatter!(
            d.df.Nthreads,
            d.df.Total,
            yerror = d.df.TotalStd,
            label = d.funcname,
            marker = m,
        )
        display(p)
        value_at_1 = filter(row -> row.Nthreads == 1, d.df).Total[1]
        color = p.subplots[1].series_list[end].plotattributes[:linecolor]
        p = plot!(1:152, x -> value_at_1 / x, label = nothing, linecolor = color) #= Adjust =#
        display(p)
        m, s = iterate(markers, s)
    end
    p = plot!(
        yaxis = :log10,
        xaxis = :log10,
        minorticks = true,
        minorgrid = true,
        leg = :outerright,
        size = (1000, 500),
    )
    ylims!(2.0e-5, 2.0e4) # Adjust
    title!("Total time by function (N=25,Nlen=14)")
    xlabel!("Nthreads")
    ylabel!("Time (s)")
    println("Saving $savename")
    savefig(savename)

    return p
end

"""
Compute the difference between the total time measured
and the sum of all the other measured times.
"""
function get_total_minus_sum(dfs)
    partial_dfs =
        [select(t.df, "Nthreads", "Total") for t in dfs if t.funcname != "Total Time"]
    total_times =
        [select(t.df, "Nthreads", "Total") for t in dfs if t.funcname == "Total Time"][1]


    # Joining all the dfs based on Nthreads
    joined = copy(partial_dfs[1][:, ["Nthreads"]])

    for d in partial_dfs
        joined = innerjoin(joined, d, on = "Nthreads", makeunique = true)
    end

    time_sums = DataFrame(
        :Nthreads => joined.Nthreads,
        :TotalSum => sum([joined[:, n] for n in names(joined, r"Total")]),
    )

    two_totals = innerjoin(total_times, time_sums, on = "Nthreads")

    DataFrame(
        :Nthreads => two_totals.Nthreads,
        :Difference => two_totals.Total - two_totals.TotalSum,
    )

end

"""
The directory should contain /only/
the data files
"""
function main(directory)
    filenames = [joinpath(directory, f) for f in readdir(directory)]
    dfs = make_dfs(filenames)
    print(
        "Difference between the total time measured\n and the sum of all the other measured times:\n",
    )
    dfs |> get_total_minus_sum |> print
    println()
    plot_dfs(dfs, "scaling-by-function.pdf")
end
