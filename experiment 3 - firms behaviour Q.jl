# Experiment 3

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex3_v1_total_surplus = []
ex3_v1_producer_surplus = []
ex3_v1_consumer_surplus = []
ex3_v1_price = []
ex3_v1_quantity_produced = []
ex3_v1_quantity_sold = []
ex3_v1_quantity_leased = []
ex3_v1_producer_surplus_singleton = []

ex3_v2_total_surplus = []
ex3_v2_producer_surplus = []
ex3_v2_consumer_surplus = []
ex3_v2_price = []
ex3_v2_quantity_produced = []
ex3_v2_quantity_sold = []
ex3_v2_quantity_leased = []
ex3_v2_producer_surplus_singleton = []

for i in 1:100

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex3_v1_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, true, 0)
    push!(ex3_v1_total_surplus, calculate_surplus(ex3_v1_sim, "total", true))
    push!(ex3_v1_producer_surplus, calculate_surplus(ex3_v1_sim, "producer", true))
    push!(ex3_v1_consumer_surplus, calculate_surplus(ex3_v1_sim, "consumer,total",true))
    push!(ex3_v1_price, calculate_price_history.(ex3_v1_sim.sellers))
    push!(ex3_v1_quantity_produced, getfield.(ex3_v1_sim.sellers, :quantity_produced_history))
    push!(ex3_v1_quantity_sold, getfield.(ex3_v1_sim.sellers, :quantity_sold_history))
    push!(ex3_v1_quantity_leased, getfield.(ex3_v1_sim.sellers, :quantity_leased_history))
    push!(ex3_v1_producer_surplus_singleton, calculate_profit_history.(ex3_v1_sim.sellers))

    ex3_v2_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, false, 0)
    push!(ex3_v2_total_surplus, calculate_surplus(ex3_v2_sim, "total", true))
    push!(ex3_v2_producer_surplus, calculate_surplus(ex3_v2_sim, "producer", true))
    push!(ex3_v2_consumer_surplus, calculate_surplus(ex3_v2_sim, "consumer,total",true))
    push!(ex3_v2_price, calculate_price_history.(ex3_v2_sim.sellers))
    push!(ex3_v2_quantity_produced, getfield.(ex3_v2_sim.sellers, :quantity_produced_history))
    push!(ex3_v2_quantity_sold, getfield.(ex3_v2_sim.sellers, :quantity_sold_history))
    push!(ex3_v2_quantity_leased, getfield.(ex3_v2_sim.sellers, :quantity_leased_history))
    push!(ex3_v2_producer_surplus_singleton, calculate_profit_history.(ex3_v2_sim.sellers))

end


ex2_p1 = plot_ecdf(ex3_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex3_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

ex2_p3 = plot_ecdf(ex3_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(ex3_v1_consumer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_consumer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex3_v1_quantity_sold)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex3_v1_quantity_leased)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_leased)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(getindex.(ex3_v1_producer_surplus_singleton,1))), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(getindex.(ex3_v2_producer_surplus_singleton,1))), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

L1 = getindex.(ex3_v12_r, 1) 
L2 = getindex.(ex3_v12_r, 2)

R1 = getindex.(ex3_v12_i, 1)
R2 = getindex.(ex3_v12_i, 2)


hm_eq = [mean_c(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[(L1 .== l1) .& (R1 .== r1)]) for l1 in sort(unique(L1)), r1 in sort(unique(R1))]'
hm_dq = [mean_c(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[(L1 .== l1) .& (R1 .== r1)]) for l1 in sort(unique(L1)), r1 in sort(unique(R1))]'

ex4_p1 = heatmap(sort(unique(L1)), sort(unique(R1)), hm_eq, xlabel = "Decrease rate", ylabel = "Increase rate", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8)
ex4_p1 = heatmap(sort(unique(L1)), sort(unique(R1)), hm_dq, xlabel = "Decrease rate", ylabel = "Increase rate", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8)

####

hm_eq = [mean_c(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[(R1 .== r1) .& (R2 .== r2)]) for r1 in sort(unique(R1)), r2 in sort(unique(R2))]'
hm_dq = [mean_c(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[(R1 .== r1) .& (R2 .== r2)]) for r1 in sort(unique(R1)), r2 in sort(unique(R2))]'

max_plot = maximum([maximum(hm_eq), maximum(hm_dq)])
min_plot = minimum([minimum(hm_eq), minimum(hm_dq)])

ex4_p1 = heatmap(sort(unique(R1)), sort(unique(R2)), hm_eq[2:end,2:end], xlabel = "Decrease rate, my", ylabel = "Decrease rate, competitors", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8, clim = (min_plot, max_plot))

ex4_p1 = heatmap(sort(unique(R1)), sort(unique(R2)), hm_dq[2:end,2:end], xlabel = "Decrease rate, my", ylabel = "Decrease rate, competitors", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość początkowa", titlefontsize = 8, clim = (min_plot, max_plot))

####


dR = R1 .- R2 # increase
dR = round.(dR, digits=2)
plot(sort(unique(dR)), [mean(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[dR .== dR_u]) for dR_u in sort(unique(dR))])
plot!(sort(unique(dR)), [mean(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[dR .== dR_u]) for dR_u in sort(unique(dR))])

dL = L1 .- L2 # decrease
dL = round.(dL, digits=2)
plot(sort(unique(dL)), [mean(mean.(getindex.(ex3_v1_producer_surplus_singleton, 1))[dL .== dL_u]) for dL_u in sort(unique(dL))])
plot!(sort(unique(dL)), [mean(mean.(getindex.(ex3_v2_producer_surplus_singleton, 1))[dL .== dL_u]) for dL_u in sort(unique(dL))])


ex3_p1 = plot_ecdf(ex3_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex3_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

ex3_p2 = plot_ecdf(ex3_v1_consumer_surplus, "equal Q", "Consumer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_consumer_surplus, "not equal Q", "Consumer Surplus", "Probability", "ECDF", false)

ex3_p3 = plot_ecdf(ex3_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex3_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex3_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_sold)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)



mean.(mean.(ex3_v1_price))

ex3_p3 = plot_ecdf(mean.(mean.(ex3_v1_price)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)