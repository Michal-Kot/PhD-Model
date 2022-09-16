# Experiment 1

include(pwd() * "\\methods\\methods.jl")

# Impact on social communication of simulation dynamics

ex1_v1_total_surplus = []
ex1_v1_producer_surplus = []
ex1_v1_consumer_surplus = []
ex1_v1_price = []
ex1_v1_quantity_produced = []
ex1_v1_quantity_sold = []
ex1_v1_quantity_leased = []
ex1_v1_producer_surplus_singleton = []

ex1_v2_total_surplus = []
ex1_v2_producer_surplus = []
ex1_v2_consumer_surplus = []
ex1_v2_price = []
ex1_v2_quantity_produced = []
ex1_v2_quantity_sold = []
ex1_v2_quantity_leased = []
ex1_v2_producer_surplus_singleton = []

ex1_v12_λ = []

for i in 1:1000

    if (mod(i,10) == 0) | (i == 1)
        println(i)
    end

    λ_ind = rand()
    λ_wom = rand()

    push!(ex1_v12_λ, [λ_ind, λ_wom])

    ex1_v1_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", λ_ind, λ_wom, "stochastic", [0.1, 0.1], [1.1, 1.1], [[0.5, 1.5], [0.5, 1.5]], [[0.2, 0.6], [0.2, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, 0)

    push!(ex1_v1_total_surplus, calculate_surplus(ex1_v1_sim, "total", true))
    push!(ex1_v1_producer_surplus, calculate_surplus(ex1_v1_sim, "producer", true))
    push!(ex1_v1_consumer_surplus, calculate_surplus(ex1_v1_sim, "consumer,total",true))
    push!(ex1_v1_price, calculate_price_history.(ex1_v1_sim.sellers))
    push!(ex1_v1_quantity_produced, getfield.(ex1_v1_sim.sellers, :quantity_produced_history))
    push!(ex1_v1_quantity_sold, getfield.(ex1_v1_sim.sellers, :quantity_sold_history))
    push!(ex1_v1_quantity_leased, getfield.(ex1_v1_sim.sellers, :quantity_leased_history))
    push!(ex1_v1_producer_surplus_singleton, calculate_profit_history.(ex1_v1_sim.sellers))

    ex1_v2_sim = TO_GO(250, 2, 200, 300, [1.0, 1.0], [0.5, 0.5], [1.0, 1.0], [0.3, 0.3], "random", λ_ind, λ_wom, "stochastic",[0.1, 0.1], [1.1, 1.1], [[0.9, 1.5], [0.6, 1.1]], [[0.2, 0.5], [0.3, 0.6]],[[0.8,2.], [0.8,2.]], 0.1, true, 0)
    push!(ex1_v2_total_surplus, calculate_surplus(ex1_v2_sim, "total", true))
    push!(ex1_v2_producer_surplus, calculate_surplus(ex1_v2_sim, "producer", true))
    push!(ex1_v2_consumer_surplus, calculate_surplus(ex1_v2_sim, "consumer,total",true))
    push!(ex1_v2_price, calculate_price_history.(ex1_v2_sim.sellers))
    push!(ex1_v2_quantity_produced, getfield.(ex1_v2_sim.sellers, :quantity_produced_history))
    push!(ex1_v2_quantity_sold, getfield.(ex1_v2_sim.sellers, :quantity_sold_history))
    push!(ex1_v2_quantity_leased, getfield.(ex1_v2_sim.sellers, :quantity_leased_history))
    push!(ex1_v2_producer_surplus_singleton, calculate_profit_history.(ex1_v2_sim.sellers))

end


L1 = getindex.(ex1_v12_λ, 1) 
L2 = getindex.(ex1_v12_λ, 2)

L1, L1u = cut_integer(L1, 5)
L2, L2u = cut_integer(L2, 5)

#### Total market surplus

hm_eq = [mean_c(ex1_v1_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'
hm_dq = [mean_c(ex1_v2_total_surplus[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

max_plot = maximum([maximum(hm_eq), maximum(hm_dq)])
min_plot = minimum([minimum(hm_eq), minimum(hm_dq)])

ex4_p1 = heatmap(L1u, L2u, hm_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, identyczna jakość", titlefontsize = 8)

savefig(ex4_p1, pwd() * "\\plots\\ex1_equal quality heatmap lambdas.svg")

ex4_p2 = heatmap(L1u, L2u, hm_dq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka całkowita, różna jakość", titlefontsize = 8)

savefig(ex4_p2, pwd() * "\\plots\\ex1_different quality heatmap lambdas.svg")

#%% Producer 1 surplus, better producer

hm_better_eq = [mean_c(sum.(getindex.(ex1_v1_producer_surplus_singleton, 1))[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'
hm_better_dq = [mean_c(sum.(getindex.(ex1_v2_producer_surplus_singleton, 1))[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

hm_worse_eq = [mean_c(sum.(getindex.(ex1_v1_producer_surplus_singleton, 2))[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'
hm_worse_dq = [mean_c(sum.(getindex.(ex1_v2_producer_surplus_singleton, 2))[(L1 .== l1) .& (L2 .== l2)]) for l1 in sort(unique(L1)), l2 in sort(unique(L2))]'

max_plot = maximum([maximum(hm_better_eq), maximum(hm_better_dq), maximum(hm_worse_eq), maximum(hm_worse_dq)])
min_plot = minimum([minimum(hm_better_eq), minimum(hm_better_dq),minimum(hm_worse_eq), minimum(hm_worse_dq)])

ex4_p1 = heatmap(L1u, L2u, hm_better_eq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta , identyczna jakość", titlefontsize = 7, clim = (min_plot, max_plot))

savefig(ex4_p1, pwd() * "\\plots\\ex1_equal quality heatmap lambdas b eq.svg")

ex4_p2 = heatmap(L1u, L2u, hm_better_dq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta, wyższa jakość", titlefontsize = 8, clim = (min_plot, max_plot))

savefig(ex4_p2, pwd() * "\\plots\\ex1_different quality heatmap lambdas b dq.svg")

ex4_p2 = heatmap(L1u, L2u, hm_worse_dq, xlabel = "λ_ind, produkty oceniane osobiście", ylabel = "λ_wom, produkty oceniane przez sąsiadów", title = "Współczynniki siły nowych sygnałów a nadwyżka producenta, niższa jakość", titlefontsize = 8, clim = (min_plot, max_plot))

savefig(ex4_p2, pwd() * "\\plots\\ex1_different quality heatmap lambdas w dq.svg")