include(pwd() * "\\methods\\methods.jl")

#################################### AUX FUNCTIONS ##############################################################

sim_single = TO_GO(250, 2, 300, 500, [0.4, 0.4], [1.25, 1.25], "random", 0.25, 0.25, "stochastic", 1.1, [[0.25, 0.75], [0.25, 0.75]], [[0.2, 0.8], [0.2, 0.8]], [[0.0, 2.], [0.0, 2.]], 0.5, true, true, 1)

StatsPlots.groupedbar(hcat(getfield.(sim_single.sellers, :reselling_history)...), bar_position = :stack, linecolor = nothing)

Plots.plot(getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected,1) .== 1) .& (getindex.(sim_single.profit_expected,2) .== "p")],3))
Plots.plot!(calculate_profit_history(sim_single.sellers[1])[3:end])

Plots.plot(getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected,1) .== 2) .& (getindex.(sim_single.profit_expected,2) .== "p")],3))
Plots.plot!(calculate_profit_history(sim_single.sellers[2])[3:end])

Plots.plot(getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected,1) .== 1) .& (getindex.(sim_single.profit_expected,2) .== "ed")],3))
Plots.plot!(getindex.(sim_single.profit_expected[(getindex.(sim_single.profit_expected,1) .== 1) .& (getindex.(sim_single.profit_expected,2) .== "qp")],3))
Plots.plot!(sim_single.sellers[1].quantity_produced_history)
Plots.plot!(sim_single.sellers[1].quantity_sold_history .+ sim_single.sellers[1].quantity_leased_history)

Plots.plot(calculate_profit_history(sim_single.sellers[1])[3:end])

Plots.plot(calculate_price_history.(sim_single.sellers))

ex4_p1 = Plots.plot(calculate_profit_history.(sim_single.sellers), color = ["blue"  "orange"], xlabel = "t", ylabel = "Nadwyżka producenta", label = ["Producent 1" "Producent 2"], title = "Nadwyżka producenta")

Plots.savefig(ex4_p1, pwd() * "\\plots\\ex4_prod surplus.svg")

ex4_p21 = plot_quantity(sim_single.sellers,1)
ex4_p22 = plot_quantity(sim_single.sellers,2)

StatsPlots.groupedbar([sim_single.sellers[1].quantity_sold_history .+ sim_single.sellers[1].quantity_leased_history sim_single.sellers[2].quantity_sold_history .+ sim_single.sellers[2].quantity_leased_history], bar_position = :stack, linecolor = nothing)

Plots.plot(sum(getfield.(sim_single.sellers, :quantity_produced_history)))

Plots.savefig(ex4_p21, pwd() * "\\plots\\ex4_prod quant 1.svg")
Plots.savefig(ex4_p22, pwd() * "\\plots\\ex4_prod quant 2.svg")

Plots.plot(sim_single.sellers[1].quality_history, color = "blue", linewidth = 2, xlabel = "t", ylabel = "Jakość / oczekiwana jakość", label = "Producent 1 - jakość")
Plots.plot!(mean([getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)]), color = "blue", linewidth = 2, label = "Producent 1 - oczekiwana jakość", linestyle = :dot)
Plots.plot!(sim_single.sellers[2].quality_history, color = "orange", linewidth = 2, label = "Producent 2 - jakość")
Plots.plot!(mean([getindex.(x,2) for x in getfield.(sim_single.buyers, :quality_expectation_history)]), color = "orange", linewidth = 2, label = "Producent 2 - oczekiwana jakość", linestyle = :dot)

Plots.plot(sim_single.sellers[1].durability_history, color = "blue", linewidth = 2, xlabel = "t", ylabel = "Jakość / oczekiwana trwałość", label = "Producent 1 - trwałość")
Plots.plot!(mean([getindex.(x,1) for x in getfield.(sim_single.buyers, :durability_expectation_history)]), color = "blue", linewidth = 2, label = "Producent 1 - oczekiwana trwałość", linestyle = :dot)
Plots.plot!(sim_single.sellers[2].durability_history, color = "orange", linewidth = 2, label = "Producent 2 - trwałość")
Plots.plot!(mean([getindex.(x,2) for x in getfield.(sim_single.buyers, :durability_expectation_history)]), color = "orange", linewidth = 2, label = "Producent 2 - oczekiwana trwałość", linestyle = :dot)

Plots.plot(sim_single.sellers[1].margin_history, color = "blue", linewidth = 2, xlabel = "t", ylabel = "Marża", label = "Producent 1 - marża")
Plots.plot!(sim_single.sellers[2].margin_history, color = "orange", linewidth = 2, label = "Producent 2 - marża")


Plots.savefig(ex4_p3, pwd() * "\\plots\\ex4_qual dura exp.svg")

ex4_p4 = Plots.plot(sim_single.sellers[1].margin_history .- 1, color = "blue", xlabel = "t", ylabel = "Marża %", label = "Producent 1", title = "Marża %")
Plots.plot!(sim_single.sellers[2].margin_history .- 1, color = "orange", label = "Producent 2")

Plots.savefig(ex4_p4, pwd() * "\\plots\\ex4_margin.svg")


ex4_p5 = Plots.plot(calculate_price_history(sim_single.sellers[1]), color = "blue", xlabel = "t", ylabel = "Cena", label = "Producent 1", title = "Cena")
Plots.plot!(calculate_price_history(sim_single.sellers[2]), color = "orange", label = "Producent 2")

Plots.savefig(ex4_p5, pwd() * "\\plots\\ex4_price.svg")

Plots.plot(getindex.(mean([b.quality_expectation_history for b in sim_single.buyers]),1) .- calculate_price_history(sim_single.sellers[1]))
Plots.plot!(getindex.(mean([b.quality_expectation_history for b in sim_single.buyers]),2) .- calculate_price_history(sim_single.sellers[2]))

any([mean(x[(end-100):end]) <= 1 for x in getfield.(sim_single.sellers, :quantity_produced_history)])

ex4_v1_total_surplus = []
ex4_v1_producer_surplus = []
ex4_v1_consumer_surplus = []
ex4_v1_price = []
ex4_v1_quality = []
ex4_v1_durability = []
ex4_v1_margin = []
ex4_v1_quantity_produced = []
ex4_v1_quantity_sold = []
ex4_v1_quantity_leased = []
ex4_v1_producer_surplus_singleton = []
ex4_v1_reselling = []
ex4_v1_consumer_surplus_pm = []
ex4_v1_consumer_surplus_sm_s = []
ex4_v1_consumer_surplus_sm_b = []
ex4_v1_market_dead = []
ex4_v1_quality = []
ex4_v1_durability = []
ex4_v1_quality_exp = []
ex4_v1_durability_exp = []

ex4_v2_total_surplus = []
ex4_v2_producer_surplus = []
ex4_v2_consumer_surplus = []
ex4_v2_price = []
ex4_v2_quality = []
ex4_v2_durability = []
ex4_v2_margin = []
ex4_v2_quantity_produced = []
ex4_v2_quantity_sold = []
ex4_v2_quantity_leased = []
ex4_v2_producer_surplus_singleton = []
ex4_v2_reselling = []
ex4_v2_consumer_surplus_pm = []
ex4_v2_consumer_surplus_sm_s = []
ex4_v2_consumer_surplus_sm_b = []
ex4_v2_market_dead = []
ex4_v2_quality = []
ex4_v2_durability = []
ex4_v2_quality_exp = []
ex4_v2_durability_exp = []

ex4_v12_nn = []

for i in 1:400

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    nn = sample(LinRange(0:25:400))

    push!(ex4_v12_nn, nn)

    ex4_v1_sim = TO_GO(800, 2, 200, nn, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1., 1.], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.8,2.], [0.8,2.]], 0.1, true, true, )

    push!(ex4_v1_total_surplus, calculate_surplus(ex4_v1_sim, "total", false))
    push!(ex4_v1_producer_surplus, calculate_surplus(ex4_v1_sim, "producer", false))
    push!(ex4_v1_consumer_surplus, calculate_surplus(ex4_v1_sim, "consumer,total", false))
    push!(ex4_v1_price, calculate_price_history.(ex4_v1_sim.sellers))
    push!(ex4_v1_quality, getfield.(ex4_v1_sim.sellers, :quality_history))
    push!(ex4_v1_durability, getfield.(ex4_v1_sim.sellers, :durability_history))
    push!(ex4_v1_margin, getfield.(ex4_v1_sim.sellers, :margin_history))
    push!(ex4_v1_quantity_produced, getfield.(ex4_v1_sim.sellers, :quantity_produced_history))
    push!(ex4_v1_quantity_sold, getfield.(ex4_v1_sim.sellers, :quantity_sold_history))
    push!(ex4_v1_quantity_leased, getfield.(ex4_v1_sim.sellers, :quantity_leased_history))
    push!(ex4_v1_producer_surplus_singleton, calculate_profit_history.(ex4_v1_sim.sellers))
    push!(ex4_v1_reselling, getfield.(ex4_v1_sim.sellers, :reselling_history))
    push!(ex4_v1_consumer_surplus_pm, calculate_surplus(ex4_v1_sim, "consumer,pm", false))
    push!(ex4_v1_consumer_surplus_sm_b, calculate_surplus(ex4_v1_sim, "consumer,sm,b", false))
    push!(ex4_v1_consumer_surplus_sm_s, calculate_surplus(ex4_v1_sim, "consumer,sm,s", false))
    push!(ex4_v1_market_dead, [mean(x[(end-100):end]) <= 1 for x in getfield.(ex4_v1_sim.sellers, :quantity_produced_history)])
    push!(ex4_v1_quality, getfield.(ex4_v1_sim.sellers, :quality_history))
    push!(ex4_v1_durability, getfield.(ex4_v1_sim.sellers, :durability_history))
    push!(ex4_v1_quality_exp, mean([b.quality_expectation_history for b in ex4_v1_sim.buyers]))
    push!(ex4_v1_durability_exp, mean([b.durability_expectation_history for b in ex4_v1_sim.buyers]))

    ex4_v2_sim = TO_GO(800, 2, 200, nn, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[1.0, 1.5], [0.5, 1.0]], [[0.3, 0.6], [0.2, 0.5]], [[.8, 2.], [.8, 2.]], 0.10, true, true, 0)
    push!(ex4_v2_total_surplus, calculate_surplus(ex4_v2_sim, "total", false))
    push!(ex4_v2_producer_surplus, calculate_surplus(ex4_v2_sim, "producer", false))
    push!(ex4_v2_consumer_surplus, calculate_surplus(ex4_v2_sim, "consumer,total",false))
    push!(ex4_v2_price, calculate_price_history.(ex4_v2_sim.sellers))
    push!(ex4_v2_quality, getfield.(ex4_v2_sim.sellers, :quality_history))
    push!(ex4_v2_durability, getfield.(ex4_v2_sim.sellers, :durability_history))
    push!(ex4_v2_margin, getfield.(ex4_v2_sim.sellers, :margin_history))
    push!(ex4_v2_quantity_produced, getfield.(ex4_v2_sim.sellers, :quantity_produced_history))
    push!(ex4_v2_quantity_sold, getfield.(ex4_v2_sim.sellers, :quantity_sold_history))
    push!(ex4_v2_quantity_leased, getfield.(ex4_v2_sim.sellers, :quantity_leased_history))
    push!(ex4_v2_producer_surplus_singleton, calculate_profit_history.(ex4_v2_sim.sellers))
    push!(ex4_v2_reselling, getfield.(ex4_v2_sim.sellers, :reselling_history))
    push!(ex4_v2_consumer_surplus_pm, calculate_surplus(ex4_v2_sim, "consumer,pm", false))
    push!(ex4_v2_consumer_surplus_sm_b, calculate_surplus(ex4_v2_sim, "consumer,sm,b", false))
    push!(ex4_v2_consumer_surplus_sm_s, calculate_surplus(ex4_v2_sim, "consumer,sm,s", false))
    push!(ex4_v2_market_dead, [mean(x[(end-100):end]) <= 1 for x in getfield.(ex4_v2_sim.sellers, :quantity_produced_history)])
    push!(ex4_v2_quality, getfield.(ex4_v2_sim.sellers, :quality_history))
    push!(ex4_v2_durability, getfield.(ex4_v2_sim.sellers, :durability_history))
    push!(ex4_v2_quality_exp, mean([b.quality_expectation_history for b in ex4_v2_sim.buyers]))
    push!(ex4_v2_durability_exp, mean([b.durability_expectation_history for b in ex4_v2_sim.buyers]))

end

Plots.plot(sort(unique(ex4_v12_nn)), [mean(all.(ex4_v1_market_dead)[ex4_v12_nn .== x]) for x in sort(unique(ex4_v12_nn))])




Plots.plot(sim_single.sellers[1].margin_history)

Plots.plot()

Plots.plot(sum_of_geom_series.(getfield.(sim_single.sellers, :quality_history)[1], getfield.(sim_single.sellers, :durability_history)[1]), linewidth = 2, xlabel = "T", ylabel = "Użyteczność / oczekiwana jakość", label = "Użyteczność, producent 1", legend = :outerbottom)
Plots.plot!(sum_of_geom_series.(getindex.(mean([b.quality_expectation_history for b in sim_single.buyers]), 1), getindex.(mean([b.durability_expectation_history for b in sim_single.buyers]), 1)) , label = "Oczekiwana użyteczność, producent 1")
Plots.plot!(sum_of_geom_series.(getfield.(sim_single.sellers, :quality_history)[2], getfield.(sim_single.sellers, :durability_history)[2]), linewidth = 2, label = "Użyteczność, producent 2")
Plots.plot!(sum_of_geom_series.(getindex.(mean([b.quality_expectation_history for b in sim_single.buyers]), 2), getindex.(mean([b.durability_expectation_history for b in sim_single.buyers]), 2)), label = "Oczekiwana użyteczność, producent 2")

plot_phasediagram(sim_single.sellers[1])

Plots.scatter( sim_single.sellers[1].quality_history .- mean([getindex.(x,1) for x in getfield.(sim_single.buyers, :quality_expectation_history)]), sim_single.sellers[1].quantity_produced_history, alpha = collect(1:1000) ./ 1000, markerstrokewidth = 0)

Plots.scatter( sim_single.sellers[2].quality_history .- mean([getindex.(x,2) for x in getfield.(sim_single.buyers, :quality_expectation_history)]), sim_single.sellers[2].quantity_produced_history)

plot_quantity(sim_single)

Plots.scatter(sum_of_geom_series.(getfield.(sim_single.sellers, :quality_history)[1], getfield.(sim_single.sellers, :durability_history)[1]) .- sum_of_geom_series.(getindex.(mean([b.quality_expectation_history for b in sim_single.buyers]), 1), getindex.(mean([b.durability_expectation_history for b in sim_single.buyers]), 1)), calculate_profit_history.(sim_single.sellers)[1], smooth = true)

Plots.plot(cumsum.(calculate_profit_history.(sim_single.sellers)))
sum([count(getfield.(x, :d) .== "s") for x in getfield.(sim_single.buyers, :unit_buying_selling_history)])

sum(sum(getfield.(sim_single.sellers, :reselling_history)))


Plots.plot(getfield.(sim_single.sellers, :reselling_history))

calculate_surplus(sim_single, "producer", true)

plot(getfield.(sim_single.sellers, :selling_income_history))
plot!(getfield.(sim_single.sellers, :leasing_income_history))

plot(getfield.(sim_single.sellers, :expected_income_history))
plot!(getfield.(sim_single.sellers, :selling_income_history) .+ getfield.(sim_single.sellers, :leasing_income_history))

sum.(getfield.(sim_single.sellers, :selling_income_history) .+ getfield.(sim_single.sellers, :leasing_income_history).- 1 .* getfield.(sim_single.sellers, :cost_of_production_history))

plot(getfield.(sim_single.sellers, :selling_income_history) .+ getfield.(sim_single.sellers, :leasing_income_history))

plot!(-1 .* getfield.(sim_single.sellers, :cost_of_production_history))

plot(getfield.(sim_single.sellers, :selling_income_history) .+ getfield.(sim_single.sellers, :leasing_income_history) .- 1 .* getfield.(sim_single.sellers, :cost_of_production_history))

plot(getfield.(sim_single.sellers, :selling_income_history) .+ getfield.(sim_single.sellers, :leasing_income_history) -1 .* getfield.(sim_single.sellers, :cost_of_production_history))


plot(getfield.(sim_single.sellers, :durability_history))
plot!(getfield.(sim_single.sellers, :quality_history))
plot(getfield.(sim_single.sellers, :margin_history))

plot(getfield.(sim_single.sellers, :quantity_produced_history))
plot!(getfield.(sim_single.sellers, :quantity_leased_history))
plot!(getfield.(sim_single.sellers, :quantity_sold_history))

plot(getfield.(sim_single.sellers, :quantity_produced_history) .- getfield.(sim_single.sellers, :quantity_leased_history) .- getfield.(sim_single.sellers, :quantity_sold_history))

plot(calculate_profit_history.(sim_single.sellers))

plot_quantity(sim_single)

plot(calculate_profit_history.(sim_single.sellers))

scatter(getfield.(sim_single.sellers, :margin_history), calculate_profit_history.(sim_single.sellers), smooth = true)

plot(getfield.(sim_single.sellers, :utilization_cost_history))

sum(sim_single.buyers[2].realized_surplus_pm_history)
sum(sim_single.buyers[7].realized_surplus_sm_s_history)
sum(sim_single.buyers[2].realized_surplus_sm_b_history)

sim_single.sellers[1].utlization_cost_history

plot(calculate_profit_history.(sim_single.sellers))
plot(calculate_price_history.(sim_single.sellers))

calculate_cost_history.(sim_single.sellers)

plot_quantity(sim_single.sellers)



plot(calculate_price_history.(sim_single.sellers))
plot(calculate_profit_history.(sim_single.sellers))


plot(calculate_surplus(sim_single, "total", false))
plot!(calculate_surplus(sim_single, "producer", false))
plot!(calculate_surplus(sim_single, "consumer,pm", false))
plot!(calculate_surplus(sim_single, "consumer,sm,b", false))
plot!(calculate_surplus(sim_single, "consumer,sm,s", false))

plot(getfield.(sim_single.sellers, :reselling_history))

sum.([any.(x) for x in getfield.(buyers, :unit_possessed_history)])

plot(getindex.(getindex(getfield.(buyers, :surplus_history), 1), 1))





plot_quantity(sim_single.sellers, 1)

length(sim_single.sellers)

#######################################

s = collect(LinRange(0:0.01:1.0))

k = [1.0, 0.8, 1.2]
d = [0.8, 0.5, 0.9]
c = [0.4, 0.3, 0.5]
m = [1.0, 1.0, 1.0]

ke = k ./ (1 .- d)
p = c .* m .* ke

sur = [s .* kei .- pi for (kei, pi) in zip(ke, p)]

plot(s,sur, ylim = (0, Inf))

choice(x) = any(x .>= 0) ? argmax(x) : 0



plot!(s, choice.([getindex.(sur, x) for x in 1:length(s)]), xlabel = "Skłonność do zapłaty za jakość")

