# Experiment 4

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex4_v1_total_surplus = []
ex4_v1_producer_surplus = []
ex4_v1_consumer_surplus = []
ex4_v1_price = []
ex4_v1_quantity_produced = []
ex4_v1_quantity_sold = []
ex4_v1_producer_surplus_singleton = []

ex4_v2_total_surplus = []
ex4_v2_producer_surplus = []
ex4_v2_consumer_surplus = []
ex4_v2_price = []
ex4_v2_quantity_produced = []
ex4_v2_quantity_sold = []
ex4_v2_producer_surplus_singleton = []

ex4_v3_total_surplus = []
ex4_v3_producer_surplus = []
ex4_v3_consumer_surplus = []
ex4_v3_price = []
ex4_v3_quantity_produced = []
ex4_v3_quantity_sold = []
ex4_v3_producer_surplus_singleton = []

ex4_v4_total_surplus = []
ex4_v4_producer_surplus = []
ex4_v4_consumer_surplus = []
ex4_v4_price = []
ex4_v4_quantity_produced = []
ex4_v4_quantity_sold = []
ex4_v4_producer_surplus_singleton = []

for i in 1:500

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex4_v1_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [0.1, 0.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.,2.], [0.,2.]], 0.50, true)
    push!(ex4_v1_total_surplus, calculate_surplus(ex4_v1_sim, "total", true))
    push!(ex4_v1_producer_surplus, calculate_surplus(ex4_v1_sim, "producer", true))
    push!(ex4_v1_consumer_surplus, calculate_surplus(ex4_v1_sim, "consumer,total",true))
    push!(ex4_v1_price, calculate_price_history.(ex4_v1_sim.sellers))
    push!(ex4_v1_quantity_produced, getfield.(ex4_v1_sim.sellers, :quantity_produced_history))
    push!(ex4_v1_quantity_sold, getfield.(ex4_v1_sim.sellers, :quantity_sold_history))
    push!(ex4_v1_producer_surplus_singleton, calculate_profit_history.(ex4_v1_sim.sellers))

    ex4_v2_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [0.1, 0.1], [[1.0, 1.5], [0.5, 1.0]], [[0.4, 0.6], [0.2, 0.4]], [[0.,2.], [0.,2.]], 0.50, true)
    push!(ex4_v2_total_surplus, calculate_surplus(ex4_v2_sim, "total", true))
    push!(ex4_v2_producer_surplus, calculate_surplus(ex4_v2_sim, "producer", true))
    push!(ex4_v2_consumer_surplus, calculate_surplus(ex4_v2_sim, "consumer,total",true))
    push!(ex4_v2_price, calculate_price_history.(ex4_v2_sim.sellers))
    push!(ex4_v2_quantity_produced, getfield.(ex4_v2_sim.sellers, :quantity_produced_history))
    push!(ex4_v2_quantity_sold, getfield.(ex4_v2_sim.sellers, :quantity_sold_history))
    push!(ex4_v2_producer_surplus_singleton, calculate_profit_history.(ex4_v2_sim.sellers))

    ex4_v3_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [0.1, 0.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]], [[0.,2.], [0.,2.]], 0.50, false)
    push!(ex4_v3_total_surplus, calculate_surplus(ex4_v3_sim, "total", true))
    push!(ex4_v3_producer_surplus, calculate_surplus(ex4_v3_sim, "producer", true))
    push!(ex4_v3_consumer_surplus, calculate_surplus(ex4_v3_sim, "consumer,total",true))
    push!(ex4_v3_price, calculate_price_history.(ex4_v3_sim.sellers))
    push!(ex4_v3_quantity_produced, getfield.(ex4_v3_sim.sellers, :quantity_produced_history))
    push!(ex4_v3_quantity_sold, getfield.(ex4_v3_sim.sellers, :quantity_sold_history))
    push!(ex4_v3_producer_surplus_singleton, calculate_profit_history.(ex4_v3_sim.sellers))


    ex4_v4_sim = TO_GO(150, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.4, 0.4], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [0.1, 0.1], [[1.0, 1.5], [0.5, 1.0]], [[0.4, 0.6], [0.2, 0.4]], [[0.,2.], [0.,2.]], 0.50, false)
    push!(ex4_v4_total_surplus, calculate_surplus(ex4_v4_sim, "total", true))
    push!(ex4_v4_producer_surplus, calculate_surplus(ex4_v4_sim, "producer", true))
    push!(ex4_v4_consumer_surplus, calculate_surplus(ex4_v4_sim, "consumer,total",true))
    push!(ex4_v4_price, calculate_price_history.(ex4_v4_sim.sellers))
    push!(ex4_v4_quantity_produced, getfield.(ex4_v4_sim.sellers, :quantity_produced_history))
    push!(ex4_v4_quantity_sold, getfield.(ex4_v4_sim.sellers, :quantity_sold_history))
    push!(ex4_v4_producer_surplus_singleton, calculate_profit_history.(ex4_v4_sim.sellers))


end

ex4_p1 = plot_ecdf(ex4_v1_total_surplus, "Identyczna jakość, rynek wtórny istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex4_v2_total_surplus, "Różna jakość, rynek wtórny istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)
plot_ecdf(ex4_v3_total_surplus, "Identyczna jakość, rynek wtórny nie istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)
plot_ecdf(ex4_v4_total_surplus, "Różna jakość, rynek wtórny nie istnieje", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)


ex4_p2 = plot_ecdf(ex4_v1_consumer_surplus, "Identyczna jakość, rynek wtórny istnieje", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex4_v2_consumer_surplus, "Różna jakość, rynek wtórny istnieje", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex4_v3_consumer_surplus, "Identyczna jakość, rynek wtórny nie istnieje", "Consumer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex4_v4_consumer_surplus, "Różna jakość, rynek wtórny nie istnieje", "Consumer Surplus", "Probability", "ECDF", false)

ex4_p3 = plot_ecdf(ex4_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex4_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex4_v3_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(ex4_v4_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex4_p3 = plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v1_producer_surplus_singleton, 1)]), 2.5), "Identyczna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", true)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v2_producer_surplus_singleton, 1)]),2.5), "Różna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v3_producer_surplus_singleton, 1)]),2.5), "Identyczna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v4_producer_surplus_singleton, 1)]),2.5), "Różna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o wyższej jakości", "Probability", "ECDF", false)

ex4_p3 = plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v1_producer_surplus_singleton, 1)]), 2.5), "Identyczna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", true)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v2_producer_surplus_singleton, 2)]),2.5), "Różna jakość, rynek wtórny istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v3_producer_surplus_singleton, 2)]),2.5), "Identyczna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", false)
plot_ecdf(trim_outliers(mean.([x[51:150] for x in getindex.(ex4_v4_producer_surplus_singleton, 2)]),2.5), "Różna jakość, rynek wtórny nie istnieje", "Nadwyżka producenta dobra o niższej jakości", "Probability", "ECDF", false)

sec_mar_better = mean.([x[51:150] for x in getindex.(ex4_v2_producer_surplus_singleton, 1)]) .- mean.([x[51:150] for x in getindex.(ex4_v4_producer_surplus_singleton, 1)])
sec_mar_worse = mean.([x[51:150] for x in getindex.(ex4_v2_producer_surplus_singleton, 2)]) .- mean.([x[51:150] for x in getindex.(ex4_v4_producer_surplus_singleton, 2)])

UnequalVarianceTTest(sec_mar_better, sec_mar_worse)

ex4_p3 = plot_ecdf(sum.(sum.(ex4_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex4_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex4_p3 = plot_ecdf(sum.(sum.(ex4_v1_quantity_sold)) / 200, "Identyczna jakość, rynek wtórny istnieje", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex4_v2_quantity_sold)) / 200, "Różna jakość, rynek wtórny istnieje", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(sum.(sum.(ex4_v3_quantity_sold)) / 200, "Identyczna jakość, rynek wtórny nie istnieje", "Producer Surplus", "Probability", "ECDF", false)
plot_ecdf(sum.(sum.(ex4_v4_quantity_sold)) / 200, "Różna jakość, rynek wtórny nie istnieje", "Producer Surplus", "Probability", "ECDF", false)





UnequalVarianceTTest(mean.(getindex.(ex4_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex4_v2_producer_surplus_singleton, 1)))

UnequalVarianceTTest(mean.(getindex.(ex4_v1_producer_surplus_singleton, 1)), mean.(getindex.(ex4_v3_producer_surplus_singleton, 1)))

UnequalVarianceTTest(mean.(getindex.(ex4_v2_producer_surplus_singleton, 1)), mean.(getindex.(ex4_v4_producer_surplus_singleton, 1)))

mean.(getindex.(ex4_v2_producer_surplus_singleton, 1))
mean.(getindex.(ex4_v3_producer_surplus_singleton, 1))
mean.(getindex.(ex4_v4_producer_surplus_singleton, 1))
mean.(getindex.(ex4_v4_producer_surplus_singleton, 1))

mean.(mean.(ex4_v1_price))

ex4_p3 = plot_ecdf(mean.(mean.(ex4_v1_price)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex4_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)