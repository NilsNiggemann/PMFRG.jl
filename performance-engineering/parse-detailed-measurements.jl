
import Statistics: mean, std
import DataFrames: DataFrame, select, filter
using Plots


function parse_file(filename)
    open(filename) do file
        lines = readlines(file)

        detailed_measurement_lines = [l for l in lines if contains(l, "tag:64")]

        timed_func_names = Set([split(l)[1] for l in detailed_measurement_lines])

        measurements = [
            (func_name = func_name,
             measurements = [parse(Float32, split(l)[3])
                             for l in detailed_measurement_lines
                             if split(l)[1] == func_name])
            for func_name in timed_func_names]

        measurement_stats = [
            (funcname = row.func_name,
             mean = mean(row.measurements),
             std = std(row.measurements))
            for row in measurements]

        return (stats = measurement_stats,
                func_names = timed_func_names)
    end
end

get_nthreads(filename) =  parse(Int32,match(r"detailed-timings?-(\d+)\.out", filename)[1])

function parse_list(filenames)

    all_parses = [ (fname = fname, stats_and_funcnames = parse_file(fname))
                   for fname in filenames ]

    return ( fnames = all_parses[1].stats_and_funcnames.func_names,
             rows =
             [(Nthreads = get_nthreads(parse.fname),
               FuncName = row.funcname,
               Mean = row.mean,
               Std = row.std)
               for parse in all_parses
               for row in parse.stats_and_funcnames.stats ])
end


function make_dfs(filenames)

    (;fnames, rows) = parse_list(filenames)

	df = DataFrame(rows)

	dfs = [ (funcname = funcname,
             df = select(filter(row -> row["FuncName"] == funcname,
                                df),
                         "Nthreads", "Mean", "Std"))
            for funcname in fnames ]
end

function plot_dfs(dfs, savename)

    plot()

	for d in dfs
	    p = scatter!(d.df.Nthreads,
                     d.df.Mean,
                     yerror = d.df.Std,
                     label = d.funcname)
	    display(p)
        value_at_1 = filter(row -> row.Nthreads == 1, d.df).Mean[1]
        color = p.subplots[1].series_list[end].plotattributes[:linecolor]
        p = plot!( 1:152 #= Adjust =#, x -> value_at_1 / x, label = nothing, linecolor = color )
	    display(p)
 	end
	p = plot!(yaxis = :log10, xaxis=:log10, minorticks = true, minorgrid= true)
	ylims!(1.0e-4,10) # Adjust
    savefig(savename)
    return p
end

"""
The directory should contain /only/
the data files
"""
function main(directory)
    filenames = [ joinpath(directory,f) for f in readdir(directory) ]
    dfs = make_dfs(filenames)
    plot_dfs(dfs,"scaling-by-function.png")
end
