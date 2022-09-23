include(pwd() * "\\methods\\methods.jl")

#################################### AUX FUNCTIONS ##############################################################

sim_single = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", 0.25, 0.25, "stochastic",[0.1, 0.1], [1.1, 1.1], [[0.9, 1.5], [0.6, 1.1]], [[0.2, 0.5], [0.3, 0.6]],[[1.0,2.], [1.0,2.]], 0.0, true, true, 0)

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

