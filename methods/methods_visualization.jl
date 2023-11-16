#######################################################################################################

# Auxiliary functions used to visualize simulation results

#######################################################################################################

function plot_ecdf(is_new, metric, label; legend = :bottomright, xlabel="", ylabel="", title="", xlim=nothing, ylim=nothing, legendfontsize = 8)

    """
    Prints ECDF chart
    Inputs:
        is_new - if chart is new or shall be plotted as additional series
        metric - chosen metric to be plotted
        label - series label
    Returns:
    ECDF plot
    """

    if is_new
        Plots.plot(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=legend, xlim=xlim, ylim=ylim, legendfontsize = legendfontsize)
    else
        Plots.plot!(sort(metric), (1:length(metric))./length(metric), label = label)
    end
end
