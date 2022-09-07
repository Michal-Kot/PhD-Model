#### Visuals

function plot_quantity(sim_res)
    groupedbar(reduce(hcat,[getfield.(sim_res.sellers, :quantity_history)[x] for x in 1:length(sim_res.sellers)])[2:end,:],linecolor=nothing, bar_position = :stack, xlabel = "Time", ylabel = "Quantity of sales", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_quality_expectation(sim_res)
    plot([getindex.(mean(getfield.(sim_res.buyers, :quality_expectation_history)),x) for x in 1:sim_res.function_args.num_sellers], label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_profit_history(sim_res)
    plot(calculate_profit_history.(sim_res.sellers))
end

function plot_advertising(sim_res)
    plot(getfield.(sim_res.sellers, :advertising_history), xlabel = "Time", ylabel = "Advertising intensity [share of population targeted]", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
end

function plot_margin(sim_res, res_perc)
    if res_perc
        Plots.plot([getfield.(sim_res.sellers, :margin_history)[x] for x in 1:length(sim_res.sellers)], xlabel = "Time", ylabel = "Margin [%]", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
    else
        prices = calculate_price_history.(sim_res.sellers)
        mcs = getfield.(sim_res.sellers, :average_cost) .* getfield.(sim_res.sellers, :quality) .* getfield.(sim_res.sellers, :durability)
        margin = [p .- mc for (p, mc) in zip(prices, mcs)]
        Plots.plot(margin, xlabel = "Time", ylabel = "Margin [abs]", label = reshape("Producer " .* string.(1:sim_res.function_args.num_sellers), 1, sim_res.function_args.num_sellers))
    end
end

function plot_ecdf(metric, label, xlabel, ylabel, title, is_new)
    if is_new
        Plots.plot(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=:bottomright, legendfontsize = 6)
    else
        Plots.plot!(sort(metric), (1:length(metric))./length(metric), xlabel = xlabel, ylabel = ylabel, title = title, label = label, legend=:bottomright)
    end
end

function add_smoothing_spline(x,y,clr,λ=0.05)
    spl = fit(SmoothingSpline, x, y, λ)
    y_hat = predict(spl)
    plot!(sort(x), y_hat[sortperm(x)], label = "", color = clr)
end

function plot_phasediagram(_seller::seller, function_args::Vector)
    x = calculate_price_history(_seller)
    y = calculate_profit_history(_seller, function_args.γ, function_args.num_buyers)
    d_x = diff(x)
    push!(d_x, 0)
    d_y = diff(y)
    push!(d_y,0)
    labs = string.(1:length(x))
    points_labels = ([mod(x,100) == 0 ? labs[x] : "" for x in 1:length(labs)])
    points_labels[1] = "1"
    quiver(x,y,quiver=(d_x, d_y), color = :green)
    scatter!(x,y, xlabel = "Price", ylabel = "Profit", markershape = :none, markercolor = :white, markerstrokecolor = :white, series_annotations = points_labels, markersize = 0)
end

function plot_quantity(sellers, random_start = 50)

    quantity_produced = getfield.(sim_single.sellers, :quantity_produced_history)
    quantity_produced = [q[random_start:end] for q in quantity_produced]

    quantity_sold = getfield.(sim_single.sellers, :quantity_sold_history)
    quantity_sold = [q[random_start:end] for q in quantity_sold]

    plot(quantity_produced)
    plot!(quantity_sold)

end
