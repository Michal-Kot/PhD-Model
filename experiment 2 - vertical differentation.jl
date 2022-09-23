# Experiment 2

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex2_v1_total_surplus = []
ex2_v1_producer_surplus = []
ex2_v1_consumer_surplus = []
ex2_v1_price = []
ex2_v1_quantity_produced = []
ex2_v1_quantity_sold = []
ex2_v1_quantity_leased = []
ex2_v1_producer_surplus_singleton = []

ex2_v2_total_surplus = []
ex2_v2_producer_surplus = []
ex2_v2_consumer_surplus = []
ex2_v2_price = []
ex2_v2_quantity_produced = []
ex2_v2_quantity_sold = []
ex2_v2_quantity_leased = []
ex2_v2_producer_surplus_singleton = []

ex2_v3_total_surplus = []
ex2_v3_producer_surplus = []
ex2_v3_consumer_surplus = []
ex2_v3_price = []
ex2_v3_quantity_produced = []
ex2_v3_quantity_sold = []
ex2_v3_quantity_leased = []
ex2_v3_producer_surplus_singleton = []

ex2_v4_total_surplus = []
ex2_v4_producer_surplus = []
ex2_v4_consumer_surplus = []
ex2_v4_price = []
ex2_v4_quantity_produced = []
ex2_v4_quantity_sold = []
ex2_v4_quantity_leased = []
ex2_v4_producer_surplus_singleton = []

for i in 1:2000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    nn = sample(1:500)

    ex2_v1_sim = TO_GO(250, 2, 200, nn, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.8,2.], [0.8,2.]], 0.0, true, true, 0)
    push!(ex2_v1_total_surplus, calculate_surplus(ex2_v1_sim, "total", true))
    push!(ex2_v1_producer_surplus, calculate_surplus(ex2_v1_sim, "producer", true))
    push!(ex2_v1_consumer_surplus, calculate_surplus(ex2_v1_sim, "consumer,total",true))
    push!(ex2_v1_price, calculate_price_history.(ex2_v1_sim.sellers))
    push!(ex2_v1_quantity_produced, getfield.(ex2_v1_sim.sellers, :quantity_produced_history))
    push!(ex2_v1_quantity_sold, getfield.(ex2_v1_sim.sellers, :quantity_sold_history))
    push!(ex2_v1_quantity_leased, getfield.(ex2_v1_sim.sellers, :quantity_leased_history))    
    push!(ex2_v1_producer_surplus_singleton, calculate_profit_history.(ex2_v1_sim.sellers))


    ex2_v2_sim = TO_GO(250, 2, 200, nn, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.8,2.], [0.8,2.]], 0.0, true, false, 0)
    push!(ex2_v2_total_surplus, calculate_surplus(ex2_v2_sim, "total", true))
    push!(ex2_v2_producer_surplus, calculate_surplus(ex2_v2_sim, "producer", true))
    push!(ex2_v2_consumer_surplus, calculate_surplus(ex2_v2_sim, "consumer,total",true))
    push!(ex2_v2_price, calculate_price_history.(ex2_v2_sim.sellers))
    push!(ex2_v2_quantity_produced, getfield.(ex2_v2_sim.sellers, :quantity_produced_history))
    push!(ex2_v2_quantity_sold, getfield.(ex2_v2_sim.sellers, :quantity_sold_history))
    push!(ex2_v2_quantity_leased, getfield.(ex2_v2_sim.sellers, :quantity_leased_history))    
    push!(ex2_v2_producer_surplus_singleton, calculate_profit_history.(ex2_v2_sim.sellers))

    ex2_v3_sim = TO_GO(250, 2, 200, nn, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.8,2.], [0.8,2.]], 0.0, false, true, 0)
    push!(ex2_v3_total_surplus, calculate_surplus(ex2_v3_sim, "total", true))
    push!(ex2_v3_producer_surplus, calculate_surplus(ex2_v3_sim, "producer", true))
    push!(ex2_v3_consumer_surplus, calculate_surplus(ex2_v3_sim, "consumer,total",true))
    push!(ex2_v3_price, calculate_price_history.(ex2_v3_sim.sellers))
    push!(ex2_v3_quantity_produced, getfield.(ex2_v3_sim.sellers, :quantity_produced_history))
    push!(ex2_v3_quantity_sold, getfield.(ex2_v3_sim.sellers, :quantity_sold_history))
    push!(ex2_v3_quantity_leased, getfield.(ex2_v3_sim.sellers, :quantity_leased_history))    
    push!(ex2_v3_producer_surplus_singleton, calculate_profit_history.(ex2_v3_sim.sellers))

    ex2_v4_sim = TO_GO(250, 2, 200, nn, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.8,2.], [0.8,2.]], 0.0, false, false, 0)
    push!(ex2_v4_total_surplus, calculate_surplus(ex2_v4_sim, "total", true))
    push!(ex2_v4_producer_surplus, calculate_surplus(ex2_v4_sim, "producer", true))
    push!(ex2_v4_consumer_surplus, calculate_surplus(ex2_v4_sim, "consumer,total",true))
    push!(ex2_v4_price, calculate_price_history.(ex2_v4_sim.sellers))
    push!(ex2_v4_quantity_produced, getfield.(ex2_v4_sim.sellers, :quantity_produced_history))
    push!(ex2_v4_quantity_sold, getfield.(ex2_v4_sim.sellers, :quantity_sold_history))
    push!(ex2_v4_quantity_leased, getfield.(ex2_v4_sim.sellers, :quantity_leased_history))    
    push!(ex2_v4_producer_surplus_singleton, calculate_profit_history.(ex2_v4_sim.sellers))

end

ex2_p1 = plot_ecdf(ex2_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex2_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

ex2_p3 = plot_ecdf(ex2_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex2_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(sum.(sum.(ex2_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex2_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex2_v1_quantity_sold)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex2_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)



mean.(mean.(ex2_v1_price))

ex2_p3 = plot_ecdf(mean.(mean.(ex2_v1_price)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex2_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)