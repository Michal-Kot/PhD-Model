include(pwd() * "\\methods\\methods.jl")

#################################### AUX FUNCTIONS ##############################################################

sim_single = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", 0.25, 0.25, "stochastic", ["stochastic", "stochastic"], [[0.5, 1.0], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.,1.1], [0.,1.1]])

plot(calculate_profit_history.(sim_single.sellers))
plot(calculate_price_history.(sim_single.sellers))

calculate_cost_history.(sim_single.sellers)

plot(getfield.(sim_single.sellers, :quantity_produced_history))
plot!(getfield.(sim_single.sellers, :quantity_sold_history))



sim_single.sellers[2].durability_history

plot(calculate_price_history.(sim_single.sellers))
plot(calculate_profit_history.(sim_single.sellers))


calculate_surplus(sim_single, "producer", true)

plot(calculate_surplus(sim_single, "producer", false))
plot!(calculate_surplus(sim_single, "consumer,pm", false))
plot!(calculate_surplus(sim_single, "consumer,sm,b", false))
plot!(calculate_surplus(sim_single, "consumer,sm,s", false))

plot(getfield.(sim_single.sellers, :reselling_history))

sum.([any.(x) for x in getfield.(buyers, :unit_possessed_history)])

plot(getindex.(getindex(getfield.(buyers, :surplus_history), 1), 1))






plot(getfield.(sim_single.sellers, :durability_history))
plot(getfield.(sim_single.sellers, :quality_history))
plot(getfield.(sim_single.sellers, :margin_history))
