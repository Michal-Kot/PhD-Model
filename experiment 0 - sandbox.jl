include(pwd() * "\\methods\\methods.jl")

#################################### AUX FUNCTIONS ##############################################################

sim_single = TO_GO(500, 3, 200, 300, [1.0, 1.0, 1.0], [0.5, 0.5, 0.5], [1.0, 1.0, 1.0], [0.3, 0.3, 0.3], "random", 0.25, 0.25, "stochastic", [0.05, 0.1, 0.1], [0.15, 0.15, 0.15], [[1.0, 1.5], [0.5, 1.0], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6], [0.2, 0.6]], [[0.5,1.5], [0.5,1.5], [0.5, 1.5]], 0.25, true)

buyers = sim_single.buyers

num_buyers = 200
num_sellers = 3





sim_single.buyers[1].unit_possessed_history

sum(sim_single.buyers[2].realized_surplus_pm_history)
sum(sim_single.buyers[7].realized_surplus_sm_s_history)
sum(sim_single.buyers[2].realized_surplus_sm_b_history)

sim_single.sellers[1].utlization_cost_history

plot(calculate_profit_history.(sim_single.sellers))
plot(calculate_price_history.(sim_single.sellers))

calculate_cost_history.(sim_single.sellers)

plot_quantity(sim_single.sellers)

plot_quantity(sim_single.sellers)



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
plot(calculate_profit_history.(sim_single.sellers))

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

