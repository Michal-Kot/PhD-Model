# Experiment 3

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex3_v1_total_surplus = []
ex3_v1_producer_surplus = []
ex3_v1_consumer_surplus = []
ex3_v1_price = []
ex3_v1_quality = []
ex3_v1_durability = []
ex3_v1_margin = []
ex3_v1_quantity_produced = []
ex3_v1_quantity_sold = []
ex3_v1_quantity_leased = []
ex3_v1_producer_surplus_singleton = []
ex3_v1_reselling = []
ex3_v1_consumer_surplus_pm = []
ex3_v1_consumer_surplus_sm_s = []
ex3_v1_consumer_surplus_sm_b = []

ex3_v2_total_surplus = []
ex3_v2_producer_surplus = []
ex3_v2_consumer_surplus = []
ex3_v2_price = []
ex3_v2_quality = []
ex3_v2_durability = []
ex3_v2_margin = []
ex3_v2_quantity_produced = []
ex3_v2_quantity_sold = []
ex3_v2_quantity_leased = []
ex3_v2_producer_surplus_singleton = []
ex3_v2_reselling = []
ex3_v2_consumer_surplus_pm = []
ex3_v2_consumer_surplus_sm_s = []
ex3_v2_consumer_surplus_sm_b = []

for i in 1:2000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    ex3_v1_sim = TO_GO(250, 2, 200, 300, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, true, 0)

    push!(ex3_v1_total_surplus, calculate_surplus(ex3_v1_sim, "total", false))
    push!(ex3_v1_producer_surplus, calculate_surplus(ex3_v1_sim, "producer", false))
    push!(ex3_v1_consumer_surplus, calculate_surplus(ex3_v1_sim, "consumer,total", false))
    push!(ex3_v1_price, calculate_price_history.(ex3_v1_sim.sellers))
    push!(ex3_v1_quality, getfield.(ex3_v1_sim.sellers, :quality_history))
    push!(ex3_v1_durability, getfield.(ex3_v1_sim.sellers, :durability_history))
    push!(ex3_v1_margin, getfield.(ex3_v1_sim.sellers, :margin_history))
    push!(ex3_v1_quantity_produced, getfield.(ex3_v1_sim.sellers, :quantity_produced_history))
    push!(ex3_v1_quantity_sold, getfield.(ex3_v1_sim.sellers, :quantity_sold_history))
    push!(ex3_v1_quantity_leased, getfield.(ex3_v1_sim.sellers, :quantity_leased_history))
    push!(ex3_v1_producer_surplus_singleton, calculate_profit_history.(ex3_v1_sim.sellers))
    push!(ex3_v1_reselling, getfield.(ex3_v1_sim.sellers, :reselling_history))
    push!(ex3_v1_consumer_surplus_pm, calculate_surplus(ex3_v1_sim, "consumer,pm", false))
    push!(ex3_v1_consumer_surplus_sm_b, calculate_surplus(ex3_v1_sim, "consumer,sm,b", false))
    push!(ex3_v1_consumer_surplus_sm_s, calculate_surplus(ex3_v1_sim, "consumer,sm,s", false))

    ex3_v2_sim = TO_GO(250, 2, 200, 300, [0.5, 0.5], [1.0, 1.0], "random", 0.25, 0.25,  "stochastic", [0.1, 0.1], [1.1, 1.1], [[1.0, 1.5], [0.5, 1.0]], [[0.3, 0.6], [0.2, 0.5]], [[.8, 2.], [.8, 2.]], 0.1, true, true, 0)
    push!(ex3_v2_total_surplus, calculate_surplus(ex3_v2_sim, "total", false))
    push!(ex3_v2_producer_surplus, calculate_surplus(ex3_v2_sim, "producer", false))
    push!(ex3_v2_consumer_surplus, calculate_surplus(ex3_v2_sim, "consumer,total",false))
    push!(ex3_v2_price, calculate_price_history.(ex3_v2_sim.sellers))
    push!(ex3_v2_quality, getfield.(ex3_v2_sim.sellers, :quality_history))
    push!(ex3_v2_durability, getfield.(ex3_v2_sim.sellers, :durability_history))
    push!(ex3_v2_margin, getfield.(ex3_v2_sim.sellers, :margin_history))
    push!(ex3_v2_quantity_produced, getfield.(ex3_v2_sim.sellers, :quantity_produced_history))
    push!(ex3_v2_quantity_sold, getfield.(ex3_v2_sim.sellers, :quantity_sold_history))
    push!(ex3_v2_quantity_leased, getfield.(ex3_v2_sim.sellers, :quantity_leased_history))
    push!(ex3_v2_producer_surplus_singleton, calculate_profit_history.(ex3_v2_sim.sellers))
    push!(ex3_v2_reselling, getfield.(ex3_v2_sim.sellers, :reselling_history))
    push!(ex3_v2_consumer_surplus_pm, calculate_surplus(ex3_v2_sim, "consumer,pm", false))
    push!(ex3_v2_consumer_surplus_sm_b, calculate_surplus(ex3_v2_sim, "consumer,sm,b", false))
    push!(ex3_v2_consumer_surplus_sm_s, calculate_surplus(ex3_v2_sim, "consumer,sm,s", false))

end

Plots.scatter(mean.(getindex.(ex3_v1_margin, 1)) .- 1, mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), markeralpha = 0.25, markerstrokewidth = 0, label = "Identyczna jakość", xlabel = "Marżowość", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v1_margin, 1)) .- 1, mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), "blue")

Plots.scatter(mean.(getindex.(ex3_v2_margin, 1)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), markeralpha = 0.25, markerstrokewidth = 0, label = "Różna jakość, producent 1", xlabel = "Marża %", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 1)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), "blue")

Plots.scatter(mean.(getindex.(ex3_v2_margin, 2)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), markeralpha = 0.25, markerstrokewidth = 0, label = "Różna jakość, producent 2", xlabel = "Marża %", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 2)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), "blue")

Plots.plot(xlabel = "Marża %", ylabel = "Zysk")
add_smoothing_spline(mean.(getindex.(ex3_v1_margin, 1)) .- 1, mean.(getindex.(ex3_v1_producer_surplus_singleton, 1)), "red", "Identyczna jakość")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 1)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 1)), "green", "Różna jakość, producent 1")
add_smoothing_spline(mean.(getindex.(ex3_v2_margin, 2)) .- 1, mean.(getindex.(ex3_v2_producer_surplus_singleton, 2)), "blue", "Różna jakość, producent 2")

Plots.scatter([getindex.([mean.(p) for p in ex3_v2_margin],x) for x in 1:2], [getindex.([mean.(p) for p in ex3_v2_producer_surplus_singleton],x) for x in 1:2], xlim = (0.8,2.0), ylim = (-5,30))

plot([mean(sum.(getindex.(ex3_v1_producer_surplus_singleton, 1))[ex3_in .== x]) for x in sort(unique(ex3_in_d))])
plot([mean(sum.(getindex.(ex3_v2_producer_surplus_singleton, 1))[ex3_in_d .== x]) for x in sort(unique(ex3_in_d))])

scatter(getindex.(ex3_in, 1), sum.(getindex.(ex3_v1_quantity_leased, 1)))
scatter(ex3_in_d, sum.(getindex.(ex3_v2_producer_surplus_singleton, 1)))

ex2_p1 = plot_ecdf(ex3_v1_total_surplus, "Identyczna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", true)
plot_ecdf(ex3_v2_total_surplus, "Różna jakość", "Całkowita nadwyżka", "F(x)", "Dystrybuanta empiryczna - nadwyżka całkowita", false)

ex2_p3 = plot_ecdf(ex3_v1_producer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_producer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(ex3_v1_consumer_surplus, "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(ex3_v2_consumer_surplus, "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

ex2_p3 = plot_ecdf(sum.(sum.(ex3_v1_quantity_produced)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_produced)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

plot_ecdf(sum.(sum.(ex3_v1_quantity_sold)) .+ sum.(sum.(ex3_v1_quantity_leased)), "equal Q", "Producer Surplus", "Probability", "ECDF", true)
plot_ecdf(sum.(sum.(ex3_v2_quantity_sold)) .+ sum.(sum.(ex3_v2_quantity_leased)), "not equal Q", "Producer Surplus", "Probability", "ECDF", false)

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